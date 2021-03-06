/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.servlets.dazzle.datasource;

import javax.servlet.*;
import java.util.*;
import java.io.*;
import org.biojava.bio.seq.*;
import org.biojava.bio.symbol.*;
import org.biojava.utils.xml.*;

/**
 * Abstract DazzleDataSource implementation which provides default implementations
 * for many methods.  This is a useful starting point for writing many
 * data source implementations.
 *
 * @author Thomas Down
 * @version 1.00
 */

public abstract class AbstractBiojavaFeatureSource 
    extends    AbstractDazzleDataSource
    implements BiojavaFeatureSource 


{
    private ServletContext context;
    /**
     * Basic implementation of the <code>init</code> method, which just records
     * the ServletContext.  Subclasses will normally want to provide their
     * own init method, but this should be chained onto <code>super.init</code>
     * so that the <code>log</code> and <code>getServletContext</code> methods
     * work correctly.
     */
    
    public void init(ServletContext ctx)
        throws DataSourceException
    {
        this.context = ctx;
    }
    
    
    /**
     * Return the <code>ServletContext</code> which was attached when this datasource
     * was initialized.
     */
    
    public ServletContext getServletContext() {
        return context;
    }


    /**
     * Log a message via the attached ServletContext.
     */
    
    public void log(String msg) {
        getServletContext().log(msg);
    }
    
    /**
     * Log a message and exception via the attached ServletContext.
     */
    
    public void log(String msg, Throwable t) {
        getServletContext().log(msg, t);
    }
      
    /**
     * Get the specified sequence.  This is used by the AbstractDataSource <code>getFeatures</code>,
     * <code>getLandmarkLength</code>, and <code>getFeaturesByID</code>
     * methods, and is a required part of the interface for DazzleReferenceSources.
     */
    
    public abstract Sequence getSequence(String ref)
        throws DataSourceException, NoSuchElementException;
    
    /**
     * Return all the features on the requested sequence
     */

    public FeatureHolder getFeatures(String ref)
        throws DataSourceException, NoSuchElementException
    {
        return getSequence(ref);
    }

    /**
     * Return the length of the requested sequence
     */

    public int getLandmarkLength(String ref)
        throws DataSourceException, NoSuchElementException
    {
        Sequence seq = getSequence(ref);
        if (seq != null) {
            return seq.length();
        } else {
            return -1;
        }
    }

    

    /**
     * Default implementation which returns an autogenerated ID.  This
     * takes the form __dazzle__&lt;type&gt;_&lt;refseq&gt;_&lt;start&gt;_&lt;stop&gt.
     * This method should be overriden whereever possible, but it provides a useful
     * fallback for features which don't have a natural ID.
     */
    
    public String getFeatureID(Feature f) {
        StringBuffer sb = new StringBuffer();
        sb.append("__dazzle__");
        sb.append(pack(f.getType()));
        sb.append('_');
        sb.append(pack(f.getSequence().getName()));
        sb.append('_');
        sb.append(f.getLocation().getMin());
        sb.append('_');
        sb.append(f.getLocation().getMax());
        return sb.toString();
    }
    
    private String pack(String s) {
        if (s.indexOf('_') < 0) {
            return s;
        } else {
            StringBuffer sb = new StringBuffer();
            for (int i = 0; i < s.length(); ++i) {
                char c = s.charAt(i);
                if (c != '_') {
                    sb.append(c);
                } else {
                    sb.append("__");
                }
            }
            return sb.toString();
        }
    }
     
    /**
     * Default implementation which returns null
     */
     
    public String getFeatureLabel(Feature f) {
        return getFeatureID(f);
    }

    /**
     * Default implementation which returns null
     */
    
    public String getScore(Feature f) {
        return null;
    }

    /**
     * Default implementation which returns the empty map
     */
    
    public Map getLinkouts(Feature f) {
        return Collections.EMPTY_MAP;
    }

   
    /**
     * Default implementation does nothing
     */
    
    public void writeXFFDetails(XMLWriter xw, Feature f)
        throws IOException
    {
    }

    /**
     * Default implementation which returns COUNT_CALCULATE
     */

    public int countFeatures(
        String reference,
		String type
    )
	    throws DataSourceException, NoSuchElementException
    {
        return COUNT_CALCULATE;
    }

    /**
     * Default implementation which returns COUNT_CALCULATE
     */

    public int countFeatures(String reference,
			     int start,
			     int end,
			     String type)
        throws DataSourceException, NoSuchElementException
    {
        return COUNT_CALCULATE;
    }

    /**
     * Default implementation which groups by BioJava feature hierarchy.  If the
     * parent is a FEATURE instance it defines a group, otherwise no groups are
     * provided.
     */
        
    public List getGroups(Feature f) {
        FeatureHolder parent = f.getParent();
        if (parent instanceof Feature) {
            Feature fParent = (Feature) parent;
            return Collections.singletonList(
                new DASGFFGroup(
                    getFeatureID(f),
                    f.getType(),
                    getFeatureLabel(f),
                    getLinkouts(f)
                )
            );
        }
                    
        return Collections.EMPTY_LIST;
    }

    /**
     * Default implementation which returns <code>false</code>. for all features
     */
    
    public boolean getShatterFeature(Feature f) {
        return false;
    }
    
    /**
     * Default implementation which returns the ReadingFrame of FramedFeatures, otherwise <code>null</code>.
     */
     
    public FramedFeature.ReadingFrame getPhase(Feature f) {
        if (f instanceof FramedFeature) {
            return ((FramedFeature) f).getReadingFrame();
        }
        return null;
    }
    
    /**
     * Default implementation which returns the empty list
     */
    
    public List getFeatureNotes(Feature f) {
        return Collections.EMPTY_LIST;
    }
    
    /**
     * Decodes IDs generated by the defauld getFeatureID method.  Implementations
     * which provide proper features IDs should definitely override this.
     */
    
    public FeatureHolder getFeaturesByID(String id, MatchType matchType)
        throws DataSourceException
    {
        log("Decoding id " + id);
        if (id.startsWith("__dazzle__")) {
            try {
                String did = id.substring(10);
            
                String type = null;
                String seq = null;
            
                int pos = 0;
                while (type == null) {
                    pos = did.indexOf('_', pos);
                    if (did.charAt(pos + 1) != '_' || did.charAt(pos + 2) == '_') {
                        type = unpack(did.substring(0, pos));
                        pos = pos + 1;
                    } else {
                        pos = pos + 2;
                    }
                }
                int pos2 = pos;
                while (seq == null) {
                    pos2 = did.indexOf('_', pos2);
                    if (did.charAt(pos + 1) != '_' || did.charAt(pos + 2) == '_') {
                        seq = unpack(did.substring(pos, pos2));
                        pos2 = pos2 + 1;
                    } else {
                        pos2 = pos2 + 2;
                    }
                }
                int upos = did.indexOf('_', pos2);
                int start = Integer.parseInt(did.substring(pos2, upos));
                int end = Integer.parseInt(did.substring(upos + 1));
                
                log("Decoded.  type=" + type + "   seq=" + seq + "     start=" + start);
                
                FeatureHolder posFeatures = getFeatures(seq);
                return posFeatures.filter(
                    new FeatureFilter.And(
                        new FeatureFilter.ByType(type),
                        new FeatureFilter.ContainedByLocation(
                            new RangeLocation(start, end)
                        )
                    )
                );
            } catch (DataSourceException ex) {
                throw ex;
            } catch (Exception ex) {
                log("Error fetching ID " + id, ex);
                return FeatureHolder.EMPTY_FEATURE_HOLDER;
            }
        } else {
            return FeatureHolder.EMPTY_FEATURE_HOLDER;
        }
    }
    
    private String unpack(String s) {
        if (s.indexOf('_') < 0) {
            return s;
        } else {
            StringBuffer sb = new StringBuffer();
            boolean escaped = false;
            for (int i = 0; i < s.length(); ++i) {
                char c = s.charAt(i);
                if (c == '_' && !escaped) {
                    escaped = true;
                } else {
                    sb.append(c);
                    escaped = false;
                }
            }
            return sb.toString();
        }
    }
    
    /**
     * Default implementation which returns the empty set.
     */
    
    public FeatureHolder getFeaturesByGroup(String id, MatchType matchType)
        throws DataSourceException
    {
        return FeatureHolder.EMPTY_FEATURE_HOLDER;
    }
}

