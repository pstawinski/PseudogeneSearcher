package pl.genebeam.pseudogenes;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Map;

import org.asciidoctor.Asciidoctor;
import org.asciidoctor.OptionsBuilder;
import org.asciidoctor.SafeMode;

public class AdocTest {
    public static void main(String[] args) {
        String text = ".Kizmet's Favorite Authors\n" + "* Edgar Allen Poe\n" + "* Sheri S. Tepper\n" + "* Bill Bryson";

        Asciidoctor asciidoctor = Asciidoctor.Factory.create();

        Map<String, Object> options = OptionsBuilder.options().compact(false).headerFooter(true).safe(SafeMode.SAFE)
                .asMap();

        String html = asciidoctor.convert(text, options);
        System.out.println(html);
        PrintWriter printWriter;
        try {
            printWriter = new PrintWriter("/tmp/blah.html");
            printWriter.write(html);
            printWriter.close();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }

    }
}
