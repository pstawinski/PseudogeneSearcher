package pl.genebeam.utils;

public class DnaHelpers {
    private String complement(String dna) {
        StringBuilder sb = new StringBuilder(dna.length());
        for (char c : dna.toCharArray()) {
            sb.append(complement(c));
        }
        return sb.reverse().toString();
    }

    public static final char a = 'a', c = 'c', g = 'g', t = 't', n = 'n', A = 'A', C = 'C', G = 'G', T = 'T', N = 'N';

    /** Returns the complement of a single byte. */
    public char complement(final char b) {
        switch (b) {
        case a:
            return t;
        case c:
            return g;
        case g:
            return c;
        case t:
            return a;
        case A:
            return T;
        case C:
            return G;
        case G:
            return C;
        case T:
            return A;
        default:
            return b;
        }
    }
}
