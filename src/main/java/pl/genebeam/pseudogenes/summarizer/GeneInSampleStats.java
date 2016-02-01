package pl.genebeam.pseudogenes.summarizer;

import java.util.List;

import org.apache.commons.lang3.StringUtils;

import com.google.common.base.Splitter;

public class GeneInSampleStats {
    private String geneName, txName;
    private String chr;
    private int txStart, txEnd;
    private double combined, byInsertSize, byClippedSize;
    private String data;
    private String sampleName;

    private final static Splitter splitter = Splitter.on('\t').omitEmptyStrings();

    public GeneInSampleStats(String line) {
        List<String> list = splitter.splitToList(line);

        this.geneName = StringUtils.substringBefore(list.get(0), ":").intern();
        this.txName = StringUtils.substringAfter(list.get(0), ":").intern();
        this.chr = StringUtils.substringBefore(list.get(1), ":").intern();
        this.txStart = Integer.valueOf(StringUtils.substringBetween(list.get(1), ":", "-"));
        this.txEnd = Integer.valueOf(StringUtils.substringAfter(list.get(1), "-"));
        this.combined = Double.valueOf(list.get(2));
        this.byInsertSize = Double.valueOf(list.get(3));
        this.byClippedSize = Double.valueOf(list.get(4));
        this.data = list.get(5);
    }

    public double getByClippedSize() {
        return byClippedSize;
    }

    public double getByInsertSize() {
        return byInsertSize;
    }

    public double getCombined() {
        return combined;
    }

    public String getGeneName() {
        return geneName;
    }

    public String getTxName() {
        return txName;
    }

    public void setByClippedSize(double byClippedSize) {
        this.byClippedSize = byClippedSize;
    }

    public void setByInsertSize(double byInsertSize) {
        this.byInsertSize = byInsertSize;
    }

    public void setCombined(double combined) {
        this.combined = combined;
    }

    public void setGeneName(String geneName) {
        this.geneName = geneName;
    }

    public void setTxName(String txName) {
        this.txName = txName;
    }

    public String getChr() {
        return chr;
    }

    public String getData() {
        return data;
    }

    public int getTxEnd() {
        return txEnd;
    }

    public int getTxStart() {
        return txStart;
    }

    @Override
    public String toString() {
        return sampleName + ":(" + getCombined() + "," + getData() + ")";
    }

    public GeneInSampleStats sampleName(String fileName) {
        sampleName = fileName;
        return this;
    }
}
