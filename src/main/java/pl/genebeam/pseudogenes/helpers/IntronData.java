package pl.genebeam.pseudogenes.helpers;

import java.util.Optional;

import com.google.common.collect.ComparisonChain;

public class IntronData implements Comparable<IntronData> {
    private String gene, transcriptionName;
    private int index;
    private Integer readsSupportingRemovalByIS;
    private Integer readsSupportingRemovalByAlignment;

    public IntronData(int idx, int readsSupportingRemovalByIS, int readsSupportingRemovalByAlignment) {
        this.index = idx;
        this.readsSupportingRemovalByAlignment = readsSupportingRemovalByAlignment;
        this.readsSupportingRemovalByIS = readsSupportingRemovalByIS;
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
}
