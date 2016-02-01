package pl.genebeam.utils;

import com.google.common.collect.ComparisonChain;

import se.sawano.java.text.AlphanumericComparator;

public class GenomicPosition implements Comparable<GenomicPosition> {
    private final String chr;
    private final int position;

    private final static AlphanumericComparator ALPHANUMERIC_COMPARATOR = new AlphanumericComparator();

    public GenomicPosition(String chr, int position) {
        super();
        this.chr = chr;
        this.position = position;
    }

    public String getChr() {
        return chr;
    }

    public int getPosition() {
        return position;
    }

    @Override
    public int compareTo(GenomicPosition o) {
        return ComparisonChain.start().compare(chr, o.chr, ALPHANUMERIC_COMPARATOR).compare(position, o.position)
                .result();
    }

    @Override
    public String toString() {
        return chr + ":" + position;
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + ((chr == null) ? 0 : chr.hashCode());
        result = prime * result + position;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        GenomicPosition other = (GenomicPosition) obj;
        if (chr == null) {
            if (other.chr != null)
                return false;
        } else if (!chr.equals(other.chr))
            return false;
        if (position != other.position)
            return false;
        return true;
    }

}
