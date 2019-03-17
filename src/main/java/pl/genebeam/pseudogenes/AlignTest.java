package pl.genebeam.pseudogenes;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

public class AlignTest {
	private SequencePair<DNASequence, NucleotideCompound> align(String querySeq, String targetSeq) {
		try {
			DNASequence target = new DNASequence(targetSeq, AmbiguityDNACompoundSet.getDNACompoundSet());
			DNASequence query = new DNASequence(querySeq, AmbiguityDNACompoundSet.getDNACompoundSet());
			SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();

			SimpleGapPenalty gapP = new SimpleGapPenalty();
			gapP.setOpenPenalty((short) 5);
			gapP.setExtensionPenalty((short) 2);

			SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(query, target,
					PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);

			return psa;
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public static void main(String[] args) {
		AlignTest test = new AlignTest();
		SequencePair<DNASequence, NucleotideCompound> psa = test.align("AAAAATTTGGGCCC", "TTTGG");

		System.out.println(psa);
	}
}
