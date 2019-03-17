package pl.genebeam.pseudogenes.helpers;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Optional;

import org.apache.commons.collections4.CollectionUtils;

import com.fasterxml.jackson.annotation.JsonIgnore;
import com.google.common.collect.ComparisonChain;

import pl.genebeam.pseudogenes.model.AlignmentResult;

public class IntronData implements Comparable<IntronData> {
	@JsonIgnore
	private String gene, transcriptionName;
	private int index;
	private Integer readsSupportingRemovalByIS;
	private Integer readsSupportingRemovalByAlignment;
	private List<AlignmentResult> alignmentResults;

	public IntronData(int idx, int readsSupportingRemovalByIS, int readsSupportingRemovalByAlignment,
			Collection<AlignmentResult> alignmentResults) {
		this.index = idx;
		this.readsSupportingRemovalByAlignment = readsSupportingRemovalByAlignment;
		this.readsSupportingRemovalByIS = readsSupportingRemovalByIS;
		if (CollectionUtils.isNotEmpty(alignmentResults)) {
			this.alignmentResults = new ArrayList<>(alignmentResults);
		}
	}

	public String getGene() {
		return gene;
	}

	public int getIndex() {
		return index;
	}

	public Integer getReadsSupportingRemovalByAlignment() {
		return Optional.ofNullable(readsSupportingRemovalByAlignment).orElse(0);
	}

	public Integer getReadsSupportingRemovalByIS() {
		return Optional.ofNullable(readsSupportingRemovalByIS).orElse(0);
	}

	public String getTranscriptionName() {
		return transcriptionName;
	}

	public void setGene(String gene) {
		this.gene = gene;
	}

	public void setIndex(int index) {
		this.index = index;
	}

	public void setReadsSupportingRemovalByAlignment(Integer readsSupportingRemovalByAlignment) {
		this.readsSupportingRemovalByAlignment = readsSupportingRemovalByAlignment;
	}

	public void setReadsSupportingRemovalByIS(Integer readsSupportingRemovalByIS) {
		this.readsSupportingRemovalByIS = readsSupportingRemovalByIS;
	}

	public void setTranscriptionName(String transcriptionName) {
		this.transcriptionName = transcriptionName;
	}

	@Override
	public int compareTo(IntronData o) {
		return ComparisonChain.start().compare(this.gene, o.gene).compare(this.transcriptionName, o.transcriptionName)
				.compare(this.index, o.index).result();
	}

	@Override
	public String toString() {

		return "(" + index + ":" + readsSupportingRemovalByAlignment + "," + readsSupportingRemovalByIS + ")";
	}

	public List<AlignmentResult> getAlignmentResults() {
		return alignmentResults;
	}

	public void setAlignmentResults(List<AlignmentResult> alignmentResults) {
		this.alignmentResults = alignmentResults;
	}

}
