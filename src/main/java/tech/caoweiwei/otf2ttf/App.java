package tech.caoweiwei.otf2ttf;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;

import com.alibaba.fastjson.JSONObject;
import com.alibaba.fastjson.JSONReader;
import com.alibaba.fastjson.JSONWriter;

public class App {
	public static void main(String[] args) throws IOException {
		if (args.length < 1) {
			System.err.println("Need json file param!");
			System.exit(0);
		}

		System.err.println("OTF JSON file: " + args[0]);

		long beginTime = System.currentTimeMillis();

		String fileName = args[0];
		FileInputStream fis = new FileInputStream(fileName);
		JSONReader reader = new JSONReader(new BufferedReader(new InputStreamReader(fis, "UTF-8")));

		JSONWriter writer = new JSONWriter(new OutputStreamWriter(System.out, "UTF-8"));

		try {
			writer.startObject();
			reader.startObject();
			while (reader.hasNext()) {
				String key = reader.readString();
				writer.writeKey(key);

				if ("glyf".equals(key)) {
					reader.startObject();
					writer.startObject();
					while (reader.hasNext()) {
						String glyfkey = reader.readString();
						writer.writeKey(glyfkey);

						reader.startObject();
						writer.startObject();
						while (reader.hasNext()) {
							String glyfitemkey = reader.readString();
							writer.writeKey(glyfitemkey);

							if ("contours".equals(glyfitemkey)) {
								Contour[][] contours = reader.readObject(Contour[][].class);
								if (contours != null) {
									contours = new CTQ().c2qContours(contours, null, null, null);
								}
								writer.writeValue(contours);
							} else if (key.equals("stemH") || key.equals("stemV") || key.equals("hintMasks")
									|| key.equals("contourMasks")) {
								Object obj = reader.readObject(Object.class);
								writer.writeValue(null);
							} else {
								Object obj = reader.readObject(Object.class);
								writer.writeValue(obj);
							}
						}
						reader.endObject();
						writer.endObject();
					}
					reader.endObject();
					writer.endObject();
				} else if (key.equals("CFF_")) {
					reader.readObject(Object.class);
					writer.writeValue(null);
				} else if (key.equals("maxp")) {
					JSONObject obj = reader.readObject(JSONObject.class);
					obj.put("version", 1.0);
					writer.writeValue(obj);
				} else {
					Object obj = reader.readObject(Object.class);
					writer.writeValue(obj);
				}
			}
			reader.endObject();
			writer.endObject();
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			reader.close();
			writer.close();
			fis.close();
		}

		long endTime = System.currentTimeMillis();
		System.err.println("Finish [" + (endTime - beginTime) + "ms]");
	}
}
