package pl.genebeam.pseudogenes.model;

import java.util.Set;

import com.google.common.collect.ComparisonChain;

import pl.genebeam.utils.GenomicPosition;

public class PositionWithEvidence implements Comparable<PositionWithEvidence> {
	private GenomicPosition position;
	private Set<ReadSummary> reads;

	public GenomicPosition getPosition() {
		return position;
	}

	public void setPosition(GenomicPosition position) {
		this.position = position;
	}

	public Set<ReadSummary> getReads() {
		return reads;
	}

	public void setReads(Set<ReadSummary> reads) {
		this.reads = reads;
	}

	@Override
	public int compareTo(PositionWithEvidence o) {
		return ComparisonChain.start().compare(position, o.position).result();
	}

}
