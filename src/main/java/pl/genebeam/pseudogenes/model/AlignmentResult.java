package pl.genebeam.pseudogenes.model;

public class AlignmentResult {
	private String bases1;
	private String bases2;
	private int identical;
	private int length;
	private String readname;
	private String readContig;
	private int readPosition;

	public int getIdentical() {
		return identical;
	}

	public void setIdentical(int identical) {
		this.identical = identical;
	}

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	public String getBases1() {
		return bases1;
	}

	public void setBases1(String bases1) {
		this.bases1 = bases1;
	}

	public String getBases2() {
		return bases2;
	}

	public void setBases2(String bases2) {
		this.bases2 = bases2;
	}

	public String getReadname() {
		return readname;
	}

	public void setReadname(String readname) {
		this.readname = readname;
	}

	public String getReadContig() {
		return readContig;
	}

	public void setReadContig(String readContig) {
		this.readContig = readContig;
	}

	public int getReadPosition() {
		return readPosition;
	}

	public void setReadPosition(int readPosition) {
		this.readPosition = readPosition;
	}

	public AlignmentResult(String bases1, String bases2, int identical, int length, String readName, String readContig,
			int readPosition) {
		super();
		this.bases1 = bases1;
		this.bases2 = bases2;
		this.identical = identical;
		this.length = length;
		this.readname = readName;
		this.readContig = readContig;
		this.readPosition = readPosition;
	}

}
