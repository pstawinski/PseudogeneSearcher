package pl.genebeam.pseudogenes.model;

import com.fasterxml.jackson.annotation.JsonInclude;

@JsonInclude(JsonInclude.Include.NON_NULL)
public class ReadSummary {
	private String name;
	private String bases;
	private String contig;
	private Integer position;

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getBases() {
		return bases;
	}

	public void setBases(String bases) {
		this.bases = bases;
	}

	public String getContig() {
		return contig;
	}

	public void setContig(String contig) {
		this.contig = contig;
	}

	public Integer getPosition() {
		return position;
	}

	public void setPosition(Integer position) {
		this.position = position;
	}

	public ReadSummary(String name, String bases) {
		super();
		this.name = name;
		this.bases = bases;
	}

	public ReadSummary() {
		// TODO Auto-generated constructor stub
	}

	public ReadSummary(String contig, Integer position, String name, String bases) {
		super();
		this.contig = contig;
		this.position = position;
		this.name = name;
		this.bases = bases;
	}

}
