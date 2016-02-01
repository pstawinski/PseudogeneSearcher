package pl.genebeam.pseudogenes;

import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.SortedSet;
import java.util.Spliterator;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;
import java.util.zip.GZIPInputStream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Stopwatch;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Range;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import pl.genebeam.pseudogenes.helpers.Report;
import pl.genebeam.pseudogenes.helpers.Transcript;
import pl.genebeam.utils.GenomicPosition;
import pl.genebeam.utils.RangeMultimapGeneral;

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

    @Parameter(names = "--genes", description = "Genes description")
    private String genesFile;

    @Parameter(names = "--omit-duplicated-reads", description = "Omit duplicates")
    private boolean omitDuplicatedReads = false;

    @Parameter(names = "--only-soft-clipped", description = "Use data only from soft clipped reads, omit these having deletion but no soft clipped ends")
    private boolean onlySoftClipped = false;

    @Parameter(names = "--help", help = true)
    private boolean help;

    @Parameter(names = "--threads", description = "Number of threads")
    private int threads = 1;

    @Parameter(names = "--position", description = "Only selected position, in format chr9:39898200-39909240")
    private String position = null;

    @Parameter(names = "--output", description = "Ouput file, give - for stdout")
    private String output;

    private final Map<String, Transcript> transcriptsByGeneName = new HashMap<>();
    private final RangeMultimapGeneral<GenomicPosition, Transcript> transcriptsMap = new RangeMultimapGeneral<>();

    private final Multimap<Transcript, GenomicPosition> disconcordants = ArrayListMultimap.create();

    private final static GenomicPosition NOT_DEFINED = new GenomicPosition("_", 0);

    /**
     * reads with such insert size are ommited in IS analysis
     */
    private static final int TOO_BIG_INSERT_SIZE_LIMIT = 100000;

    /**
     * The cutoff where we consider insert size strange (too big)
     */
    private static final int UNUSUAL_INSERT_SIZE_LIMIT_MAX = 550;
    private static final int UNUSUAL_INSERT_SIZE_LIMIT_MIN = 120;

    private static final int READ_LENGTH = 101;

    private static final boolean SMAD4_ONLY = false;

    private static final String FIRST_DELIM = ":";
    private static final String SECOND_DELIM = "!";

    private static final int GRACE_EXON_MAPPED_TO_INTRON_BASES = 3;

    public static void main(String[] args) {
        Stopwatch stopwatch = Stopwatch.createStarted();
        log.debug("Hello World!");
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
    private SamReader bamReader;

    private SAMSequenceDictionary bamSequenceDictionary;

    private void go() throws IOException {

        referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(referenceFastaFile),
                true, true);

        bamFile = "/wum/pio/experiments/pseudogeny/enah/Sample_R_867_14_2_23_S54.marked.realigned.fixed.recal.bam";

        // position = "chr2:228,189,896-228,194,480".replace(",", "");

        loadTranscriptsData();

        SamInputResource samInputResource;

        if ("-".equals(bamFile)) {
            samInputResource = SamInputResource.of(System.in);
        } else {
            samInputResource = SamInputResource.of(new File(bamFile));
        }

        OutputStream os;
        if ("-".equals(output)) {
            os = System.out;
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
            bamStream = bamStream.parallel();
            final Stream<SAMRecord> bamStreamFinal = bamStream;
            ForkJoinPool forkJoinPool = new ForkJoinPool(threads);
            try {
                forkJoinPool.submit(() -> bamStreamFinal.forEach(this::processRead)).get();
            } catch (InterruptedException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            } catch (ExecutionException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
        } else {
            bamStream.forEach(this::processRead);
        }

        bamReader.close();

        IOUtils.write(report.toString(), os);

        StringBuilder sb = new StringBuilder();
        for (Map.Entry<Transcript, Collection<GenomicPosition>> transcriptToPositions : disconcordants.asMap()
                .entrySet()) {

            Transcript transcript = transcriptToPositions.getKey();
            sb.append(transcript.getGene() + " (" + transcript.getChr() + ":" + transcript.getTxStart() + "):\t"
                    + transcriptToPositions.getValue().stream().sorted().map(GenomicPosition::toString)
                            .collect(Collectors.joining(","))
                    + "\n");

        }
        FileUtils.write(new File("/tmp/disc.log"), sb.toString());

        os.close();

    }

    private Report report = new Report();

    private void loadTranscriptsData() throws IOException, FileNotFoundException {
        BufferedReader transcriptsReader = new BufferedReader(
                new InputStreamReader(new GZIPInputStream(new FileInputStream(genesFile))));
        try {
            transcriptsReader.lines().skip(1).forEach(this::readTranscript);
        } catch (Exception e) {
            e.printStackTrace();
        }
        transcriptsReader.close();
    }

    private int counter = 0;

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
            counter++;
            if (counter % 1_000_000 == 0) {
                String position = "";
                if (!read.getReadUnmappedFlag()) {
                    position = bamSequenceDictionary.getSequence(read.getReferenceIndex()).getSequenceName() + ":"
                            + read.getUnclippedStart();
                }
                log.debug("Counter: " + counter + "; " + position);

                // log.debug(report.toString() +
                // "\n--------------------------");
            }

            if (omitDuplicatedReads && read.getDuplicateReadFlag()) {
                return;
            }

            // if
            // ("HISEQ578:547:HBDYTADXX:1:1106:19643:23911".equals(read.getReadName()))
            // {
            // log.debug("Here");
            // }

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
                    // clipped seq analysis
                    analyzeClippedSeq(read, transcriptCoveredByRead);

                    // insert size analysis
                    analyzeInsertSize(read, transcriptCoveredByRead);

                    analyzeDisconcordant(read, transcriptCoveredByRead);
                }

            }
        } catch (Exception e) {
            log.error("Uncaught exception: ", e);
        }
    }

    private void analyzeDisconcordant(SAMRecord read, Set<Transcript> transcriptCoveredByRead) {
        if (read.getReadPairedFlag()) {
            if (!read.getMateUnmappedFlag()) {
                if (!read.getReferenceName().equals(read.getMateReferenceName())
                        || (Math.abs(read.getAlignmentStart() - read.getMateAlignmentStart()) > 20000)) {
                    GenomicPosition gp = new GenomicPosition(read.getMateReferenceName(), read.getMateAlignmentStart());
                    for (Transcript t : transcriptCoveredByRead) {
                        synchronized (disconcordants) {
                            disconcordants.put(t, gp);
                        }
                    }
                }
            } else {
                for (Transcript t : transcriptCoveredByRead) {
                    synchronized (disconcordants) {
                        disconcordants.put(t, NOT_DEFINED);
                    }
                }
            }
        }
    }

    private void analyzeClippedSeq(SAMRecord read, Collection<Transcript> transcripts) {
        // is mapped
        int unclippedStart = read.getUnclippedStart();
        int unclippedEnd = read.getUnclippedEnd();

        boolean mayContainInformationAboutIntronRemoval = false;

        // end of the soft clipped left part of the read
        Integer softClippedLeftEnd = null;
        // start of the soft clipped right part of the read
        Integer softClippedRightStart = null;

        int position = 0;
        for (CigarElement element : read.getCigar().getCigarElements()) {
            CigarOperator operator = element.getOperator();
            if (CigarOperator.SOFT_CLIP.equals(operator) || CigarOperator.DELETION.equals(operator)) {
                if (CigarOperator.SOFT_CLIP.equals(operator)) {
                    if (position == 0) {
                        softClippedLeftEnd = element.getLength();
                    } else {
                        softClippedRightStart = position;
                    }
                    mayContainInformationAboutIntronRemoval = true;
                } else if (CigarOperator.DELETION.equals(operator)) {
                    if (!onlySoftClipped)
                        mayContainInformationAboutIntronRemoval = true;
                }
            }
            position += element.getLength();
        }

        if (!mayContainInformationAboutIntronRemoval)
            return;

        for (Transcript transcript : transcripts) {
            SortedSet<Integer> coveredJunctions = transcript
                    .getJunctionsIn(Range.closedOpen(unclippedStart, unclippedEnd));

            if (!coveredJunctions.isEmpty()) {
                // read spans a junction site

                // if (read.getReadName().endsWith("1216:4199:97110"))
                // System.out.println("here");

                if (mayContainInformationAboutIntronRemoval) {
                    String readBasesString = read.getReadString();

                    for (Integer junctionReference : coveredJunctions) {
                        int junctionRead = read.getReadPositionAtReferencePosition(junctionReference, true) - 1;
                        if (junctionRead == -1) {
                            // we are out of the read, probably due to
                            // the
                            // soft clipping

                            if (transcript.isExonStartPoint(junctionReference) && softClippedLeftEnd != null) {
                                junctionRead = softClippedLeftEnd;
                            } else
                                if (transcript.isIntronStartPoint(junctionReference) && softClippedRightStart != null) {
                                junctionRead = softClippedRightStart;
                            } else {
                                // read does not support any simple
                                // theory,
                                // omit it now
                                log.trace(
                                        "According to the pseudogene theory: Shouldn't be here, we have a junction and no base at this place but no correct soft clipped part detected for read "
                                                + read.getStart() + ": " + read.getReadName());
                                continue;
                            }

                        }

                        String floatingSequence;
                        Integer intronPossiblyRemovedIndex;
                        Integer neighbourExoneIndex;
                        String neighbourExoneSequence;
                        if (transcript.isExonStartPoint(junctionReference)) {
                            floatingSequence = readBasesString.substring(0, junctionRead);
                            intronPossiblyRemovedIndex = transcript.getExonIndexAtPosition(junctionReference) - 1;
                            neighbourExoneIndex = intronPossiblyRemovedIndex;
                            if (neighbourExoneIndex == -1) {
                                // this is a first exone, possibly
                                // removed intron is "before" this exon;
                                // do omit this
                                continue;
                            }

                            Range<Integer> neighbourExoneRange = transcript.getExonWithIndex(neighbourExoneIndex);
                            neighbourExoneSequence = new String(
                                    referenceSequenceFile.getSubsequenceAt(transcript.getChr(),
                                            neighbourExoneRange.upperEndpoint() - floatingSequence.length() - 10,
                                            neighbourExoneRange.upperEndpoint()).getBases());

                        } else {
                            try {
                                floatingSequence = readBasesString.substring(junctionRead);
                            } catch (Exception e) {
                                System.out.println(e);
                                read.getReadPositionAtReferencePosition(junctionReference, true);
                                floatingSequence = "";
                            }

                            intronPossiblyRemovedIndex = transcript.getIntronIndexAtPosition(junctionReference);
                            if (intronPossiblyRemovedIndex == null) {
                                // this was the last exone, our introne
                                // is after the gene; just continue;
                                continue;
                            }
                            neighbourExoneIndex = intronPossiblyRemovedIndex + 1;
                            Range<Integer> neighbourExoneRange = transcript.getExonWithIndex(neighbourExoneIndex);

                            neighbourExoneSequence = new String(referenceSequenceFile
                                    .getSubsequenceAt(transcript.getChr(), neighbourExoneRange.lowerEndpoint(),
                                            neighbourExoneRange.lowerEndpoint() + floatingSequence.length() + 10)
                                    .getBases());
                        }

                        Range<Integer> intronPossiblyRemovedRange = transcript
                                .getIntronWithIndex(intronPossiblyRemovedIndex);
                        int intronPossiblyRemovedSize = intronPossiblyRemovedRange.upperEndpoint()
                                - intronPossiblyRemovedRange.lowerEndpoint();
                        if (floatingSequence.length() < 3 || intronPossiblyRemovedSize < 10) {
                            // omit..., too small to be considered
                        } else {

                            if (alignsToNeighourExone(floatingSequence, neighbourExoneSequence)) {

                                // let's consider
                                // intronPossiblyRemovedIndex
                                // as removed
                                log.trace("Considering intron " + intronPossiblyRemovedIndex + " as removed");

                                report.addByAlign(transcript, intronPossiblyRemovedIndex, read);
                            } else {
                                log.trace("NOT Considering intron " + intronPossiblyRemovedIndex + " as removed");
                            }

                            // log.debug(read.getUnclippedStart() + ":"
                            // +
                            // read.getReadName() + "\n" +
                            // psa.toString());
                        }
                    }
                }
            }
        }
    }

    private boolean alignsToNeighourExone(String floatingSequence, String neighbourExoneSequence) {

        if (neighbourExoneSequence.endsWith(floatingSequence) || neighbourExoneSequence.startsWith(floatingSequence)) {
            // fast check for exact alignment
            return true;
        } else {
            // do full alignment
            SequencePair<DNASequence, NucleotideCompound> psa = align(floatingSequence, neighbourExoneSequence);
            if (psa.getNumIdenticals() >= 3 && psa.getNumIdenticals() > (0.8 * floatingSequence.length())) {
                return true;
            } else {
                return false;
            }
        }
    }

    private void analyzeInsertSize(SAMRecord read, Collection<Transcript> transcripts) {
        if (read.getReadPairedFlag() && !read.getMateUnmappedFlag()) {

            // if
            // ("HISEQ:155:C61EWANXX:7:1213:15987:4728".equals(read.getReadName()))
            // {
            // System.out.println("here");
            // }

            int insertSize = Math.abs(read.getInferredInsertSize());
            String chr = bamSequenceDictionary.getSequence(read.getReferenceIndex()).getSequenceName();
            String chrMate = bamSequenceDictionary.getSequence(read.getMateReferenceIndex()).getSequenceName();

            if (insertSize != 0 && read.getMateReferenceIndex().equals(read.getReferenceIndex())) {
                if (insertSize < TOO_BIG_INSERT_SIZE_LIMIT && isUnusualInsertSize(insertSize)) {
                    Range<Integer> readRange;
                    int mateStart = read.getMateAlignmentStart();

                    if (read.getAlignmentEnd() - GRACE_EXON_MAPPED_TO_INTRON_BASES < read.getMateAlignmentStart()
                            + GRACE_EXON_MAPPED_TO_INTRON_BASES) {
                        readRange = Range.closedOpen(read.getAlignmentEnd() - GRACE_EXON_MAPPED_TO_INTRON_BASES,
                                read.getMateAlignmentStart() + GRACE_EXON_MAPPED_TO_INTRON_BASES);
                    } else if (read.getMateAlignmentStart() + READ_LENGTH
                            - GRACE_EXON_MAPPED_TO_INTRON_BASES < read.getAlignmentStart()
                                    + GRACE_EXON_MAPPED_TO_INTRON_BASES) {
                        readRange = Range.closedOpen(
                                read.getMateAlignmentStart() + READ_LENGTH - GRACE_EXON_MAPPED_TO_INTRON_BASES,
                                read.getAlignmentStart() + GRACE_EXON_MAPPED_TO_INTRON_BASES);
                    } else {
                        log.trace("Strange read starts/stops");
                        return;
                    }

                    for (Transcript transcript : transcripts) {
                        if (!transcript.getJunctionsIn(readRange).isEmpty()) {
                            // transcript.getExonesIn(readRange);
                            Collection<Integer> intronesContained = transcript.getIntronesIn(readRange);
                            int intronesSize = intronesContained.stream().map(transcript::getIntronWithIndex)
                                    .filter(Objects::nonNull).mapToInt((x) -> x.upperEndpoint() - x.lowerEndpoint())
                                    .sum();

                            // int intronesSizeReal =
                            // transcript.getIntronesSizesIn(readRange);

                            int insertSizeWithoutIntron = insertSize - intronesSize;
                            if (!isUnusualInsertSize(insertSizeWithoutIntron)) {
                                log.trace("Intrones " + intronesContained + " are considered REMOVED (" + insertSize
                                        + " > " + insertSizeWithoutIntron + ")" + read.getReadName());

                                // FIXME this
                                // intronesContained.stream().findFirst().get()
                                // is not fully correct

                                // the common case is when one of the reads does
                                // map to the part of the intron it claims to be
                                // removed... but we probably don't have enough
                                // data to know that

                                report.addByInsertSize(transcript, intronesContained.stream().findFirst().get(), read);

                            } else {
                                log.trace("Intrones " + intronesContained + " are strange, we have cannot decide ("
                                        + insertSize + " > " + insertSizeWithoutIntron + ") [" + chr + ":"
                                        + read.getStart() + "-" + read.getMateAlignmentStart() + "] "
                                        + read.getReadName());
                            }
                        }
                    }
                } else {
                    // if a pair is in the neighbour intron,
                    // consider
                    // this
                    // as a correct size and intron existing
                }
            } else {
                // some problems with mate: different chromosomes or
                // other
                // problems
                // Use this data for finding a
                if (!read.getMateReferenceIndex().equals(read.getReferenceIndex())) {

                    log.trace("SV: " + chr + ":" + read.getStart() + " vs " + chrMate + ":"
                            + read.getMateAlignmentStart());
                }
            }
        }
    }

    private SequencePair<DNASequence, NucleotideCompound> align(String querySeq, String targetSeq) {
        try {
            DNASequence target = new DNASequence(targetSeq, AmbiguityDNACompoundSet.getDNACompoundSet());
            DNASequence query = new DNASequence(querySeq, AmbiguityDNACompoundSet.getDNACompoundSet());
            SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();

            SimpleGapPenalty gapP = new SimpleGapPenalty();
            gapP.setOpenPenalty((short) 5);
            gapP.setExtensionPenalty((short) 2);

            SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(query, target,
                    PairwiseSequenceAlignerType.LOCAL, gapP, matrix);

            return psa;
        } catch (Exception e) {
            log.error("Exception", e);
        }
        return null;
    }

    private boolean isUnusualInsertSize(int insertSize) {
        return insertSize > UNUSUAL_INSERT_SIZE_LIMIT_MAX || insertSize < UNUSUAL_INSERT_SIZE_LIMIT_MIN;
    }

}
