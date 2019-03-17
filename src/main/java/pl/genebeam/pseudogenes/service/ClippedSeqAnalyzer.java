package pl.genebeam.pseudogenes.service;

import java.util.Collection;
import java.util.SortedSet;

import org.apache.commons.lang3.StringUtils;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

import com.google.common.collect.Range;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import pl.genebeam.pseudogenes.helpers.Report;
import pl.genebeam.pseudogenes.helpers.Transcript;
import pl.genebeam.pseudogenes.model.AlignmentResult;

public class ClippedSeqAnalyzer {
	private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(ClippedSeqAnalyzer.class);
	private static final boolean TRACE = true;

	public void analyzeClippedSeq(SAMRecord read, Collection<Transcript> transcripts,
			ReferenceSequenceFile referenceSequenceFile, boolean onlySoftClipped, Report report) {
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
							} else if (transcript.isIntronStartPoint(junctionReference)
									&& softClippedRightStart != null) {
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
											neighbourExoneRange.lowerEndpoint() + floatingSequence.length() + 5)
									.getBases());
						}

						Range<Integer> intronPossiblyRemovedRange = transcript
								.getIntronWithIndex(intronPossiblyRemovedIndex);
						int intronPossiblyRemovedSize = intronPossiblyRemovedRange.upperEndpoint()
								- intronPossiblyRemovedRange.lowerEndpoint();
						if (floatingSequence.length() < 5 || intronPossiblyRemovedSize < 10) {
							// omit..., too small to be considered
						} else {
							AlignmentResult ar = alignsToNeighourExone(floatingSequence, neighbourExoneSequence, read);
							if (ar != null) {

								// let's consider
								// intronPossiblyRemovedIndex
								// as removed
								if (log.isTraceEnabled()) {
									log.trace("Considering intron " + intronPossiblyRemovedIndex + " as removed, "
											+ transcript.getGene());
								}

								report.addByAlign(transcript, intronPossiblyRemovedIndex, read, ar);
							} else {
								if (log.isTraceEnabled()) {
									log.trace("NOT Considering intron " + intronPossiblyRemovedIndex + " as removed");
								}
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

	private AlignmentResult alignsToNeighourExone(String floatingSequence, String neighbourExoneSequence,
			SAMRecord read) {
		// if (neighbourExoneSequence.endsWith(floatingSequence) ||
		// neighbourExoneSequence.startsWith(floatingSequence)) {
		// // fast check for exact alignment
		// return new AlignmentResult(floatingSequence, neighbourExoneSequence,
		// floatingSequence.);
		// } else {
		// do full alignment
		SequencePair<DNASequence, NucleotideCompound> psa = align(floatingSequence, neighbourExoneSequence);
		if (psa.getNumIdenticals() >= 5 && psa.getNumIdenticals() > (0.9 * floatingSequence.length())) {
			if (log.isTraceEnabled()) {
				log.trace("Alignment: \n" + psa.toString());
			}

			String alignedString = psa.toString();

			return new AlignmentResult(StringUtils.substringBefore(alignedString, "\n"),
					StringUtils.substringAfter(alignedString, "\n").replace("\n", ""), psa.getNumIdenticals(),
					floatingSequence.length(), read.getReadName(), read.getContig(), read.getStart());
		} else {
			return null;
		}
		// }
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

}
