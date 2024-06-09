package shared;

import java.util.ArrayList;

public interface LineParser {

	/**
	 * Designate the line of text to be parsed.
	 * @param line Line to parse.
	 * @return This
	 */
	public LineParser set(byte[] line);

	/**
	 * Designate the line of text to be parsed.
	 * @param line Line to parse.
	 * @param maxTerm Stop parsing after this term.
	 * @return This
	 */
	public LineParser set(byte[] line, int maxTerm);

	/**
	 * Clear all internal state.
	 * @return This
	 */
	public LineParser clear();

	/**
	 * Reset state to the start of the current line.
	 * @return This
	 */
	public LineParser reset();

	/**
	 * Parse the designated field as an int.
	 * @param term Field to parse.
	 * @return int value.
	 */
	public int parseInt(int term);

	/**
	 * Parse the designated field as a long.
	 * @param term Field to parse.
	 * @return long value.
	 */
	public long parseLong(int term);
	
	/**
	 * Parse the designated field as a float.
	 * @param term Field to parse.
	 * @return float value.
	 */
	public float parseFloat(int term);
	
	/**
	 * Parse the designated field as a double.
	 * @param term Field to parse.
	 * @return double value.
	 */
	public double parseDouble(int term);
	
	/**
	 * Parse a byte from the designated field.
	 * @param term Field to parse.
	 * @param offset Position in the field to read.
	 * @return byte value.
	 */
	public byte parseByte(int term, int offset);
	
	/**
	 * Parse the designated field as a String.
	 * @param term Field to parse.
	 * @return A String.
	 */
	public String parseString(int term);
	
	/** Convenience method; determines whether the line starts with this term. 
	 * Implemented in Tools. */
	public boolean startsWith(String s);
	
	@Override
	public String toString();
	
	/**
	 * Format the line as a list of Strings, one per field.
	 * @return A List.
	 */
	public ArrayList<String> toList();

}