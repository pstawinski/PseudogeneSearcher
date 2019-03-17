package pl.genebeam.pseudogenes.helpers;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

public class AlignmentPairwise {
	public static void main(String[] args) throws Exception {
		String targetSeq = "CACGTTTCTTGTGGCAGCTTAAGTTTGAATGTCATTTCTTCAATGGGACGGA"
				+ "GCGGGTGCGGTTGCTGGAAAGATGCATCTATAACCAAGAGGAGTCCGTGCGCTTCGACAGC"
				+ "GACGTGGGGGAGTACCGGGCGGTGACGGAGCTGGGGCGGCCTGATGCCGAGTACTGGAACA"
				+ "GCCAGAAGGACCTCCTGGAGCAGAGGCGGGCCGCGGTGGACACCTACTGCAGACACAACTA"
				+ "CGGGGTTGGTGAGAGCTTCACAGTGCAGCGGCGAG";
		DNASequence target = new DNASequence(targetSeq, AmbiguityDNACompoundSet.getDNACompoundSet());

		String querySeq = "CGGGCCGCGGTGGACACCTACTGCAGACACAACTACGGGGTTGTGAGAGCTTCACCGTGCA" + "GCGGCGAGACGCACTCGT";
		DNASequence query = new DNASequence(querySeq, AmbiguityDNACompoundSet.getDNACompoundSet());

		SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();

		SimpleGapPenalty gapP = new SimpleGapPenalty();
		gapP.setOpenPenalty((short) 5);
		gapP.setExtensionPenalty((short) 2);

		SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(query, target,
				PairwiseSequenceAlignerType.LOCAL, gapP, matrix);

		psa.getNumIdenticals();
		System.out.println(psa);
	}
}
