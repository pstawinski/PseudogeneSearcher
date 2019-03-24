package pl.genebeam.pseudogenes;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Reader;
import java.nio.charset.StandardCharsets;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.Spliterator;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Range;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import pl.genebeam.pseudogenes.helpers.FixedBatchSpliteratorWrapper;
import pl.genebeam.pseudogenes.helpers.Report;
import pl.genebeam.pseudogenes.helpers.Transcript;
import pl.genebeam.pseudogenes.model.PositionWithEvidence;
import pl.genebeam.pseudogenes.service.ClippedSeqAnalyzer;
import pl.genebeam.pseudogenes.service.DisconcordanceAnalyzer;
import pl.genebeam.pseudogenes.service.InsertSizeAnalyzer;
import pl.genebeam.utils.GenomicPosition;
import pl.genebeam.utils.RangeMultimapGeneral;
import pl.genebeam.utils.TxNameToGeneName;

/**
 * Hello world!
 *
 */
public class App {
	private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(App.class);

	@Parameter(names = "--reference", description = "Reference fastq")
	private String referenceFastaFile;

	@Parameter(names = "--bam", description = "Input bam file")
	private String bamFile;

	@Parameter(names = "--genes", description = "Genes description from RefSeq UCSC hgTables")
	private String genesFile;

	@Parameter(names = "--pseudogenes", description = "Pseudogenes gtf description")
	private String pseudogenesFile;

	@Parameter(names = "--omit-duplicated-reads", description = "Omit duplicates")
	private boolean omitDuplicatedReads = false;

	@Parameter(names = "--only-soft-clipped", description = "Use data only from soft clipped reads, omit these having deletion but no soft clipped ends")
	private boolean onlySoftClipped = false;

	@Parameter(names = "--two-pass-run", description = "Read bam file twice, cannot be used with stdin as bam input")
	private boolean twoPass = false;

	@Parameter(names = "--help", help = true)
	private boolean help;

	@Parameter(names = "--threads", description = "Number of threads")
	private int threads = 1;

	@Parameter(names = "--position", description = "Only selected position, in format chr9:39898200-39909240")
	private String position = null;

	@Parameter(names = "--output", description = "Ouput file, give - for stdout")
	private String output;
	@Parameter(names = "--output-json", description = "Ouput file for json", required = true)
	private String jsonOutput;
	@Parameter(names = "--output-vcf", description = "Ouput file for vcf", required = true)
	private String vcfOutput;

	@Parameter(names = "--sample-name", description = "Sample name", required = false)
	private String sampleName = "name";

	@Parameter(names = "--standard-genome", description = "Value: hg19 or hg38, if given program uses it's own databases", required = false)
	private String standardGenome = null;

	private final Map<String, Transcript> transcriptsByGeneName = new HashMap<>();
	private final RangeMultimapGeneral<GenomicPosition, Transcript> transcriptsMap = new RangeMultimapGeneral<>();

	private final ClippedSeqAnalyzer clippedSeqAnalyzer = new ClippedSeqAnalyzer();
	private final InsertSizeAnalyzer insertSizeAnalyzer = new InsertSizeAnalyzer();
	private DisconcordanceAnalyzer disconcordanceAnalyzer;

	public static void main(String[] args) {
		Stopwatch stopwatch = Stopwatch.createStarted();
		log.debug("Welcome to Pseudogene Searcher!");
		App app = new App();
		JCommander cmd = new JCommander(app, args);
		try {
			app.go();
		} catch (Exception e) {
			log.error("Error:", e);
			cmd.usage();
		}
		stopwatch.stop();
		log.debug("Bye! Finished in: " + stopwatch.elapsed(TimeUnit.SECONDS) + " seconds.");

	}

	private ReferenceSequenceFile referenceSequenceFile;

	private SAMSequenceDictionary bamSequenceDictionary;
	private Report report;

	private void go() throws IOException {
		SamReader bamReader;

		referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(referenceFastaFile),
				true, true);

		Reader abysovReader = null;

		if (StringUtils.isNotBlank(standardGenome)) {
			log.info("Running with standard genome: " + standardGenome);

			switch (standardGenome.toLowerCase()) {
			case "hg19": {
				ClassLoader classLoader = getClass().getClassLoader();
				InputStream res = classLoader.getResourceAsStream("Abyzov-Genome-Res-2013-supp3.csv");
				abysovReader = new BufferedReader(new InputStreamReader(res));
				log.info("I'll load Abyzov-Genome-Res-2013-supp3 database");

				if (StringUtils.isBlank(pseudogenesFile)) {
					InputStream known = classLoader.getResourceAsStream("gencode.v19.2wayconspseudos.gtf");
					File tmp = File.createTempFile("pseudogenes", ".tmp");
					IOUtils.copy(known, new FileOutputStream(tmp));
					tmp.deleteOnExit();
					pseudogenesFile = tmp.getAbsolutePath();
					log.info("I'll load internal pseudogenes database: Gencode 19");
				}

			}
				break;
			case "hg38": {
				if (StringUtils.isBlank(pseudogenesFile)) {
					ClassLoader classLoader = getClass().getClassLoader();
					InputStream known = classLoader.getResourceAsStream("gencode.v29.2wayconspseudos.gtf");
					File tmp = File.createTempFile("pseudogenes", ".tmp");
					IOUtils.copy(known, new FileOutputStream(tmp));
					tmp.deleteOnExit();
					pseudogenesFile = tmp.getAbsolutePath();
					log.info("I'll load internal pseudogenes database: Gencode 29");
				}
			}
				break;
			default:
				throw new RuntimeException("Invalid standard genome, only hg19 and hg38 are supported");
			}

		}

		loadTranscriptsData();
		report = new Report(txNameToGeneName);

		disconcordanceAnalyzer = new DisconcordanceAnalyzer(pseudogenesFile, txNameToGeneName, abysovReader);
		IOUtils.closeQuietly(abysovReader);

		SamInputResource samInputResource;
		if ("-".equals(bamFile)) {
			samInputResource = SamInputResource.of(System.in);
		} else {
			samInputResource = SamInputResource.of(new File(bamFile));
		}

		OutputStream os;
		if ("-".equals(output)) {
			os = System.out;
			twoPass = false;
		} else {
			os = new BufferedOutputStream(new FileOutputStream(output));
		}

		bamReader = SamReaderFactory.makeDefault().open(samInputResource);

		bamSequenceDictionary = bamReader.getFileHeader().getSequenceDictionary();

		Spliterator<SAMRecord> splitIterator;
		if (StringUtils.isBlank(position)) {
			splitIterator = bamReader.spliterator();
		} else {
			log.info("Doing only a limited range analysis");
			try {
				SAMRecordIterator iterator = bamReader.queryOverlapping(StringUtils.substringBefore(position, ":"),
						Integer.valueOf(StringUtils.substringBetween(position, ":", "-")),
						Integer.valueOf(StringUtils.substringAfter(position, "-")));

				Iterable<SAMRecord> iterable = () -> iterator;
				splitIterator = iterable.spliterator();
			} catch (Exception e) {
				log.error("Unable to parse position: " + position);
				throw e;
			}
		}
		Stream<SAMRecord> bamStream = StreamSupport.stream(splitIterator, false);
		if (threads > 1) {
			bamStream = FixedBatchSpliteratorWrapper.toFixedBatchStream(bamStream.parallel(), 10_000);
		}
		bamStream.forEach(read -> processRead(read));
		bamReader.close();

		OutputStream vcfOS = new BufferedOutputStream(new FileOutputStream(vcfOutput));
		// VcfBuilder outputVcf = new VcfBuilder(vcfOS, referenceSequenceFile);

		StringBuilder sb = new StringBuilder();
		for (Map.Entry<Transcript, Collection<PositionWithEvidence>> transcriptToPositions : disconcordanceAnalyzer
				.getDisconcordants().asMap().entrySet()) {
			Transcript transcript = transcriptToPositions.getKey();

			Set<String> hitsPseudogenesOf = new HashSet<String>();
			Set<String> hitsAbyzovOf = new HashSet<String>();

			for (PositionWithEvidence pwe : transcriptToPositions.getValue()) {
				GenomicPosition gp = pwe.getPosition();
				hitsPseudogenesOf
						.addAll(disconcordanceAnalyzer.findPseudogenesIn(Range.closed(gp.move(-100), gp.move(100))));

				hitsAbyzovOf.addAll(disconcordanceAnalyzer.findHitsAbyzovOf(Range.closed(gp.move(-100), gp.move(100))));

			}

			report.addDisconcordances(transcript, hitsPseudogenesOf);
			report.addAbysovDisconcordances(transcript, hitsAbyzovOf);
			if (disconcordanceAnalyzer.isReferencePseudogene(transcript.getGene())) {
				report.addReferencePseudogene(transcript);

			}

			if (disconcordanceAnalyzer.isAbysov(transcript.getGene())) {
				report.addAbysovGene(transcript, true);

			}

			sb.append(transcript.getGene() + ":" + transcript.getName() + " (" + transcript.getChr() + ":"
					+ transcript.getTxStart() + "):\t"
					// +
					// transcriptToPositions.getValue().stream().sorted().map(PositionWithEvidence::getPosition)
					// .map(GenomicPosition::toString).collect(Collectors.joining(","))
					+ "hits: " + hitsPseudogenesOf + "\n");

			report.addPossiblePositions(transcript,
					transcriptToPositions.getValue().stream().sorted().collect(Collectors.toList()));

		}

		if (twoPass) {
			log.info("Making second pass");
			bamReader = SamReaderFactory.makeDefault().open(SamInputResource.of(new File(bamFile)));
			bamReader.forEach(read -> report.populateReadData(read));
			bamReader.close();
		}

		FileUtils.write(new File("/tmp/disc.log"), sb.toString(), StandardCharsets.UTF_8);
		IOUtils.write(report.toString(), os, StandardCharsets.UTF_8);

		// for (IdentifiedPseudogene ip : ips) {
		// outputVcf.add(ip);
		// }

		BufferedOutputStream jsonOs = new BufferedOutputStream(new FileOutputStream(jsonOutput));
		report.toVcf(referenceSequenceFile, vcfOS, jsonOs, sampleName);
		// outputVcf.close();
		os.close();

	}

	private TxNameToGeneName txNameToGeneName;

	private void loadTranscriptsData() throws IOException, FileNotFoundException {
		BufferedReader transcriptsReader = new BufferedReader(new InputStreamReader(getGenesFileStream()));
		try {
			transcriptsReader.lines().skip(1).forEach(this::readTranscript);
		} catch (Exception e) {
			e.printStackTrace();
		}
		transcriptsReader.close();

		this.txNameToGeneName = new TxNameToGeneName(getGenesFileStream());
	}

	private InputStream getGenesFileStream() throws IOException {
		InputStream genesFileStream = null;
		switch (standardGenome.toLowerCase()) {
		case "hg19": {
			ClassLoader classLoader = getClass().getClassLoader();
			log.info("I'll load Abyzov-Genome-Res-2013-supp3 database");
			if (StringUtils.isBlank(genesFile)) {
				log.info("I'll load internal hg19 genes database");
				genesFileStream = new GZIPInputStream(classLoader.getResourceAsStream("hg19.gz"));
			}
		}
			break;
		case "hg38": {
			ClassLoader classLoader = getClass().getClassLoader();
			if (StringUtils.isBlank(genesFile)) {
				log.info("I'll load internal hg38 genes database");
				genesFileStream = new GZIPInputStream(classLoader.getResourceAsStream("hg38.gz"));
			}

		}
			break;
		default:
			throw new RuntimeException("Invalid standard genome, only hg19 and hg38 are supported");
		}

		if (genesFileStream == null) {
			genesFileStream = new FileInputStream(genesFile);
		}
		return genesFileStream;
	}

	private final AtomicInteger counter = new AtomicInteger();

	private void readTranscript(String line) {
		if (StringUtils.isBlank(line)) {
			return;
		}

		Transcript transcript = new Transcript(line);
		synchronized (transcriptsByGeneName) {
			transcriptsByGeneName.put(transcript.getGene(), transcript);
			transcriptsMap.put(Range.closedOpen(new GenomicPosition(transcript.getChr(), transcript.getTxStart()),
					new GenomicPosition(transcript.getChr(), transcript.getTxEnd())), transcript);
		}

	}

	private void processRead(SAMRecord read) {
		try {
			int c = counter.incrementAndGet();
			if (c % 1_000_000 == 0) {
				String position = "";
				if (!read.getReadUnmappedFlag()) {
					position = bamSequenceDictionary.getSequence(read.getReferenceIndex()).getSequenceName() + ":"
							+ read.getUnclippedStart();
				}
				log.debug("Counter: " + c + "; " + position);
			}

			if (omitDuplicatedReads && read.getDuplicateReadFlag()) {
				return;
			}

			if (!read.getReadUnmappedFlag()) {
				String chr = bamSequenceDictionary.getSequence(read.getReferenceIndex()).getSequenceName();
				GenomicPosition readStart = new GenomicPosition(chr, read.getUnclippedStart());
				GenomicPosition readEnd = new GenomicPosition(chr, read.getUnclippedEnd());

				Set<Transcript> transcriptCoveredByRead = new HashSet<>();
				transcriptCoveredByRead.addAll(transcriptsMap.get(readStart));
				transcriptCoveredByRead.addAll(transcriptsMap.get(readEnd));

				int unclippedStart = read.getUnclippedStart();
				int unclippedEnd = read.getUnclippedEnd();
				if (unclippedStart > unclippedEnd) {
					log.warn("Strange read " + read.getReadName() + " at " + chr + ":" + unclippedStart + " with end: "
							+ chr + ":" + unclippedEnd);
					return;
				}

				if (!transcriptCoveredByRead.isEmpty()) {
					clippedSeqAnalyzer.analyzeClippedSeq(read, transcriptCoveredByRead, referenceSequenceFile,
							onlySoftClipped, report);
					insertSizeAnalyzer.analyzeInsertSize(read, transcriptCoveredByRead, bamSequenceDictionary, report);
					disconcordanceAnalyzer.analyzeDisconcordant(read, transcriptCoveredByRead);
				}
			}
		} catch (Exception e) {
			log.error("Uncaught exception: ", e);
		}
	}

}
