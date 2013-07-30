package org.biojava.servlets.dazzle.datasource;

import au.com.bytecode.opencsv.CSVReader;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.servlets.dazzle.Segment;

import javax.servlet.ServletContext;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.*;

public class TransmartAcghPlugin extends AbstractGFFFeatureSource implements DazzleReferenceSource {

    private Map seqs;
    private Set allTypes;
    private ServletContext ctx;
    private String mapMaster;
    private String fileName;

    protected List<String[]> readRegionFile(ServletContext ctx, Segment seg) throws IOException {
        BufferedReader br = new BufferedReader(new InputStreamReader(ctx.getResourceAsStream(fileName)));
        CSVReader csvReader = new CSVReader(br, '\t');
        ArrayList<String[]> results = new ArrayList<String[]>();
        try {
            String[] header;
            if ((header = csvReader.readNext()) != null) {
                results.add(header);
                String[] row;
                while ((row = csvReader.readNext()) != null) {
                    if (row.length > 3) {
                        String chromosome = row[1];
                        int start = row[2].matches("\\d+") ? Integer.valueOf(row[2]) : Integer.MIN_VALUE;
                        int stop = row[3].matches("\\d+") ? Integer.valueOf(row[3]) : Integer.MIN_VALUE;
                        if (chromosome.equals(seg.getReference())
                                && (seg.getStart() == Integer.MIN_VALUE || start >= seg.getStart())
                                && (seg.getStop() == Integer.MAX_VALUE || stop <= seg.getStop())) {
                            results.add(row);
                        }
                    }
                }
            }
        } finally {
            csvReader.close();
        }
        return results;
    }

    protected GFFFeature[] readFeatures(ServletContext ctx, Segment seg, String[] types) {
        System.out.println("Get features for " + seg + " For types: " + types);
        HashSet<String> typesSet = new HashSet<String>();
        if(types != null) {
            for (String type: types) {
                typesSet.add(type);
            }
        }
        List<String[]> rows;
        try {
            rows = readRegionFile(ctx, seg);
        } catch (IOException e) {
            e.printStackTrace();
            return new GFFFeature[0];
        }

        ArrayList<GFFFeature> features = new ArrayList<GFFFeature>();

        String[] header = rows.remove(0);
        for (String[] row : rows) {
            if(typesSet.isEmpty() || typesSet.contains("loss")) {
                GFFFeature lossFeature = new GFFFeature();
                lossFeature.setStart(row[2]);
                lossFeature.setEnd(row[3]);
                String lossUn = "loss-" + row[1] + "-" + row[2] + "-" + row[3];
                lossFeature.setLabel("lbl-" + lossUn);
                lossFeature.setName("name-" + lossUn);
                lossFeature.setType("loss");
                lossFeature.setMethod("acgh");
                lossFeature.setScore(String.valueOf((int) Math.round(Double.parseDouble(row[5]) * 100)));
                lossFeature.setPhase("-");
                lossFeature.setOrientation("0");

                features.add(lossFeature);
            }

            if(typesSet.isEmpty() || typesSet.contains("gain")) {
                GFFFeature gainFeature = new GFFFeature();
                gainFeature.setStart(row[2]);
                gainFeature.setEnd(row[3]);
                String gainUn = "loss-" + row[1] + "-" + row[2] + "-" + row[3];
                gainFeature.setLabel("lbl-" + gainUn);
                gainFeature.setName("name-" + gainUn);
                gainFeature.setType("gain");
                gainFeature.setMethod("acgh");
                gainFeature.setScore(String.valueOf((int) Math.round(Double.parseDouble(row[6]) * 100)));

                features.add(gainFeature);
            }

        }

        GFFFeature[] res = (GFFFeature[]) features.toArray(new GFFFeature[features.size()]);
        return res;
    }

    public void init(ServletContext ctx) throws DataSourceException {
        super.init(ctx);
        this.ctx = ctx;
    }

    public GFFFeature[] getFeatures(Segment seg, String[] types) throws DataSourceException {
        return readFeatures(ctx, seg, types);
    }

    public GFFFeature[] getFeatures(String reference) throws DataSourceException {
        return readFeatures(ctx, new Segment(reference), null);
    }

    public String getMapMaster() {
        return mapMaster;
    }

    public void setMapMaster(String s) {
        this.mapMaster = s;
    }

    public String getFileName() {
        return fileName;
    }

    public void setFileName(String fileName) {
        this.fileName = fileName;
    }

    /**
     * This method deals with the DAS -entry points command.
     *
     * @return a set containing the references to the entry points
     */
    public Set getEntryPoints() {
        Set<String> s = new TreeSet<String>();
        // this example has only one feature.
        // for your real data you might want to add a SQL query here.
        s.add("123");
        return s;
    }

    /**
     * This method deals with the DAS - sequence command.
     *
     * @return a biojava Sequence object
     */
    public Sequence getSequence(String ref) throws NoSuchElementException, DataSourceException {
        String seq = "ECNEUQESECNEUQESECNEUQESECNEUQESECNEUQES";

        try {
            Sequence prot = ProteinTools.createProteinSequence(seq, ref);
            return prot;
        } catch (IllegalSymbolException e) {
            throw new DataSourceException(e.getMessage());
        }
    }

}
