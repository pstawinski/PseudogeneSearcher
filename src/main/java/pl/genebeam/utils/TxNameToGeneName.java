package pl.genebeam.utils;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.nio.charset.StandardCharsets;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVRecord;

/**
 * 
 * Read gtf with ucsc_id, translate it to gene name
 * 
 * @author pio
 *
 */
public class TxNameToGeneName {

	private Map<String, String> txToGene = new HashMap<>();

	public TxNameToGeneName(InputStream stream) {
		Reader in;
		try {
			in = new InputStreamReader(stream, StandardCharsets.UTF_8.toString());

			Iterable<CSVRecord> records = CSVFormat.DEFAULT.withCommentMarker('!').withHeader().withDelimiter('\t')
					.parse(in);
			for (CSVRecord record : records) {
				String txName = record.get("#name");
				String geneName = record.get("name2");
				txToGene.put(txName, geneName);

			}
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public String getGene(String txName) {
		return txToGene.get(txName);
	}
}
