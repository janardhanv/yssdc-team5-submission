package model;

import utils.Utils;

public class GlobalContext {
	public static String gbk;
	public static int[] codes;
	public static void init(String s) {
		gbk = s;
		codes = Utils.toInt(s);
	}
}
