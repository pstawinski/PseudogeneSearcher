package pl.genebeam.pseudogenes.service;

import java.io.IOException;
import java.io.Reader;
import java.util.HashSet;
import java.util.Optional;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.lang3.StringUtils;
import org.biojava.nbio.genome.parsers.gff.FeatureI;
import org.biojava.nbio.genome.parsers.gff.FeatureList;
import org.biojava.nbio.genome.parsers.gff.GeneMarkGTFReader;
import org.biojava.nbio.genome.parsers.gff.Location;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;

import htsjdk.samtools.SAMRecord;
import pl.genebeam.pseudogenes.helpers.Transcript;
import pl.genebeam.pseudogenes.model.PositionWithEvidence;
import pl.genebeam.pseudogenes.model.ReadSummary;
import pl.genebeam.utils.GenomicPosition;
import pl.genebeam.utils.RangeMultimapGeneral;
import pl.genebeam.utils.TxNameToGeneName;

public class DisconcordanceAnalyzer {
	private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(DisconcordanceAnalyzer.class);

	private RangeMultimapGeneral<GenomicPosition, String> referencePseudogenePositions = new RangeMultimapGeneral<>();
	private RangeMultimapGeneral<GenomicPosition, String> abysovPositions = new RangeMultimapGeneral<>();
	private Set<String> referencePseudogene = new HashSet<>();
	private Set<String> abyzovGenes = new HashSet<>();
	final Multimap<Transcript, PositionWithEvidence> disconcordants = ArrayListMultimap.create();
	private final static GenomicPosition NOT_DEFINED = new GenomicPosition("_", 0);

	private final TxNameToGeneName txNameToGeneName;
	private final static Pattern COORDINATE_PATTERN = Pattern
			.compile("(?<chr>chr[0-9XYM]*):(?<start>[0-9]*)-(?<end>[0-9]*)");

	public DisconcordanceAnalyzer(String pseudogenesFile, TxNameToGeneName txNameToGeneName, Reader abysovFileReader) {
		this.txNameToGeneName = txNameToGeneName;
		if (abysovFileReader != null) {
			try {
				FeatureList listGenes = GeneMarkGTFReader.read(pseudogenesFile);
				listGenes.forEach(this::addKnownPseudogene);

				Iterable<CSVRecord> records = CSVFormat.DEFAULT.withHeader().withDelimiter('\t')
						.parse(abysovFileReader);
				for (CSVRecord record : records) {
					String category = record.get("Category");
					String geneName = record.get("Gene name");
					String coordinates = record.get("Insertion coordinate (hg19)");

					abyzovGenes.add(geneName);
					if (StringUtils.isNotBlank(coordinates)) {
						Matcher m = COORDINATE_PATTERN.matcher(coordinates);
						if (m.matches()) {
							GenomicPosition start = new GenomicPosition(m.group("chr"),
									Integer.valueOf(m.group("start")));
							GenomicPosition end = new GenomicPosition(m.group("chr"), Integer.valueOf(m.group("end")));
							abysovPositions.put(Range.closed(start, end), geneName);

						} else {
							log.warn("No match: " + coordinates);
						}
					}
				}

			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} else {
			// abysov genes will be empty
		}
	}

	private void addKnownPseudogene(FeatureI feature) {
		Location location = feature.location();
		referencePseudogenePositions.put(
				Range.closed(new GenomicPosition(feature.seqname(), Math.abs(location.getBegin())),
						new GenomicPosition(feature.seqname(), Math.abs(location.getEnd()))),
				feature.getAttribute("ucsc_id"));

		String ucscId = StringUtils.substringBefore(feature.getAttribute("ucsc_id"), ".");
		if (StringUtils.isNotBlank(txNameToGeneName.getGene(ucscId))) {
			referencePseudogene.add(txNameToGeneName.getGene(ucscId));
		}

	}

	public Set<String> findPseudogenesIn(Range<GenomicPosition> range) {
		return referencePseudogenePositions.get(range);
	}

	public Set<String> findHitsAbyzovOf(Range<GenomicPosition> range) {
		return abysovPositions.get(range);
	}

	public boolean isAbysov(String geneName) {
		return abyzovGenes.contains(geneName);
	}

	public boolean isReferencePseudogene(String geneName) {
		return referencePseudogene.contains(geneName);
	}

	public void analyzeDisconcordant(SAMRecord read, Set<Transcript> transcriptCoveredByRead) {
		if (read.getReadPairedFlag()) {
			ReadSummary readSummary = new ReadSummary();
			readSummary.setName(read.getReadName());
			readSummary.setContig(read.getContig());
			readSummary.setPosition(read.getStart());

			if (!read.getMateUnmappedFlag()) {
				if (!read.getReferenceName().equals(read.getMateReferenceName())
						|| (Math.abs(read.getAlignmentStart() - read.getMateAlignmentStart()) > 20000)) {
					GenomicPosition gp = new GenomicPosition(read.getMateReferenceName(), read.getMateAlignmentStart());
					for (Transcript t : transcriptCoveredByRead) {
						synchronized (disconcordants) {
							Optional<PositionWithEvidence> samePlace = disconcordants.get(t).stream()
									.filter(x -> x.getPosition().equals(gp)).findAny();
							if (samePlace.isPresent()) {
								samePlace.get().getReads().add(readSummary);
							} else {
								PositionWithEvidence pwe = new PositionWithEvidence();
								pwe.setPosition(gp);
								Set<ReadSummary> list = new HashSet<>();
								list.add(readSummary);
								pwe.setReads(list);
								disconcordants.put(t, pwe);
							}
						}
					}
				}
			} else {
				for (Transcript t : transcriptCoveredByRead) {
					synchronized (disconcordants) {
						// add to existing "not defined" if exists,
						// otherwise:
						// create new "not defined"
						Optional<PositionWithEvidence> notDefined = disconcordants.get(t).stream()
								.filter(x -> x.getPosition().equals(NOT_DEFINED)).findAny();
						if (notDefined.isPresent()) {

							notDefined.get().getReads().add(readSummary);
						} else {
							PositionWithEvidence pwe = new PositionWithEvidence();
							pwe.setPosition(NOT_DEFINED);
							Set<ReadSummary> list = new HashSet<>();
							list.add(readSummary);
							pwe.setReads(list);
							disconcordants.put(t, pwe);
						}
					}
				}
			}
		}
	}

	public Multimap<Transcript, PositionWithEvidence> getDisconcordants() {
		return disconcordants;
	}
}
