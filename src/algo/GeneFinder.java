package algo;

import java.util.List;

import model.CGene;

public interface GeneFinder {

	public abstract List<CGene> findGenes(String gbk);

}