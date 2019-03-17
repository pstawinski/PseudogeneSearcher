package pl.genebeam.pseudogenes.service;

import java.util.Collection;
import java.util.Objects;

import com.google.common.collect.Range;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import pl.genebeam.pseudogenes.helpers.Report;
import pl.genebeam.pseudogenes.helpers.Transcript;

public class InsertSizeAnalyzer {
	private static final org.slf4j.Logger log = org.slf4j.LoggerFactory.getLogger(InsertSizeAnalyzer.class);

	/**
	 * reads with such insert size are ommited in IS analysis
	 */
	private static final int TOO_BIG_INSERT_SIZE_LIMIT = 100000;
	private static final int GRACE_EXON_MAPPED_TO_INTRON_BASES = 3;
	private static final int READ_LENGTH = 101;
	/**
	 * The cutoff where we consider insert size strange (too big)
	 */
	private static final int UNUSUAL_INSERT_SIZE_LIMIT_MAX = 550;
	private static final int UNUSUAL_INSERT_SIZE_LIMIT_MIN = 120;

	public void analyzeInsertSize(SAMRecord read, Collection<Transcript> transcripts,
			SAMSequenceDictionary bamSequenceDictionary, Report report) {
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

	private boolean isUnusualInsertSize(int insertSize) {
		return insertSize > UNUSUAL_INSERT_SIZE_LIMIT_MAX || insertSize < UNUSUAL_INSERT_SIZE_LIMIT_MIN;
	}

}
