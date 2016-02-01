package pl.genebeam.pseudogenes.summarizer;

import java.io.IOException;
import java.nio.file.DirectoryIteratorException;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.concurrent.TimeUnit;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.google.common.base.Stopwatch;
import com.google.common.collect.ComparisonChain;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multimap;
import com.google.common.collect.Multiset;

public class SummarizeMultiple {
    private static final org.apache.log4j.Logger log = org.apache.log4j.Logger.getLogger(SummarizeMultiple.class);

    @Parameter(names = "--input", description = "Pseudogene summary files", required = false)
    private List<String> files;
    @Parameter(names = "--inputDir", description = "Pseudogene summary files", required = false)
    private String filesDir;

    @Parameter(names = "--help", help = true)
    private boolean help;

    @Parameter(names = "--output", description = "Ouput file, give - for stdout")
    private String output;

    private final static double COMMON_THRESHOLD = 0.5;
    private final static double TRUE_THRESHOLD = 0.7;

    public static void main(String[] args) {
        Stopwatch stopwatch = Stopwatch.createStarted();
        log.debug("Hello World!");
        SummarizeMultiple app = new SummarizeMultiple();
        JCommander cmd = new JCommander(app, args);
        try {
            app.go();
        } catch (Exception e) {
            log.error("Error:", e);
            cmd.usage();
        }
        stopwatch.stop();
        log.debug("Bye! Finished in: " + stopwatch.elapsed(TimeUnit.SECONDS) + " seconds.");
        System.exit(0);

    }

    private void go() {
        Path dirPath = Paths.get(filesDir);
        try (DirectoryStream<Path> stream = Files.newDirectoryStream(dirPath)) {
            for (Path entry : stream) {
                log.debug("Doing " + entry);
                Files.lines(entry).map(GeneInSampleStats::new).map(x -> x.sampleName(entry.getFileName().toString()))
                        .forEach(this::process);
            }
        } catch (DirectoryIteratorException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        geneToStats.asMap().entrySet().stream().sorted((a, b) -> ComparisonChain.start()
                .compare(commonGenes.count(a.getKey()), commonGenes.count(b.getKey())).result()).forEach(x -> {
                    System.out.println(x.getKey() + "\t" + commonGenes.count(x.getKey()) + "\t" + x.getValue());
                });
        ;

    }

    private final Multiset<String> commonGenes = HashMultiset.create();
    private final Multimap<String, GeneInSampleStats> geneToStats = HashMultimap.create();

    private void process(GeneInSampleStats stats) {
        if (stats.getCombined() >= COMMON_THRESHOLD) {
            commonGenes.add(stats.getGeneName());

            if (stats.getCombined() >= TRUE_THRESHOLD) {
                geneToStats.put(stats.getGeneName(), stats);
            }
        }

    }
}
