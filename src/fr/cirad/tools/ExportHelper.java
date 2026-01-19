package fr.cirad.tools;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Component;

@Component
public class ExportHelper {
	@Autowired private AppConfig appConfig;

	public HashMap<String, HashMap<String, String>> getOnlineOutputToolURLs() {
		HashMap<String, HashMap<String, String>> results = new LinkedHashMap<>();
		for (int i=1; ; i++)
		{
			String toolInfo = appConfig.get("onlineOutputTool_" + i);
			if (toolInfo == null)
				break;

			String[] splitToolInfo = toolInfo.split(";");
			if (splitToolInfo.length >= 2 && splitToolInfo[1].trim().length() > 0 && splitToolInfo[0].trim().length() > 0)
			{
				HashMap<String, String> aResult = new HashMap<>();
				aResult.put("url", splitToolInfo[1].trim());
				if (splitToolInfo.length >= 3 && splitToolInfo[2].trim().length() > 0)
					aResult.put("formats", splitToolInfo[2].trim());
				results.put(splitToolInfo[0].trim(), aResult);
			}
		}
		return results;
	}
	
	/**
	 * Generates customized tool URLs for a specific exported file
	 * 
	 * @param exportFormat The format of the exported file (e.g., "VCF", "PCA")
	 * @param fileUrls Map of file type to URL (e.g., "vcf" -> "http://...", "tsv" -> "http://...")
	 * @return Map with tool labels as keys and customized URLs as values
	 */
	public HashMap<String, String> getToolURLsForExport(String exportFormat, HashMap<String, String> fileUrls) {
	    HashMap<String, String> resolvedUrls = new HashMap<>();
	    HashMap<String, HashMap<String, String>> onlineTools = getOnlineOutputToolURLs();

	    for (Map.Entry<String, HashMap<String, String>> toolEntry : onlineTools.entrySet()) {
	        String toolLabel = toolEntry.getKey();
	        HashMap<String, String> toolInfo = toolEntry.getValue();

	        String urlTemplate = toolInfo.get("url");
	        String supportedFormats = toolInfo.get("formats");

	        // Check if this tool supports the export format
	        if (supportedFormats != null && !supportedFormats.isEmpty()) {
	            String[] formats = supportedFormats.split(",");
	            boolean formatSupported = false;
	            for (String format : formats) {
	                if (format.trim().equalsIgnoreCase(exportFormat)) {
	                    formatSupported = true;
	                    break;
	                }
	            }
	            if (!formatSupported) {
	                continue; // Skip this tool
	            }
	        }

	        // Track which extensions will be used by this tool
	        List<String> usedExtensions = new ArrayList<>();

	        // Replace placeholders in URL template
	        String resolvedUrl = urlTemplate;
	        Pattern placeholderPattern = Pattern.compile("\\{([^}]+)\\}");
	        Matcher matcher = placeholderPattern.matcher(urlTemplate);

	        while (matcher.find()) {
	            String placeholder = matcher.group(1);
	            String[] extensions = placeholder.split("\\|");

	            // Find matching file URL
	            String replacementUrl = null;
	            String matchedExt = null;
	            for (String ext : extensions) {
	                ext = ext.trim();
	                if (fileUrls.containsKey(ext)) {
	                    replacementUrl = fileUrls.get(ext);
	                    matchedExt = ext;
	                    break;
	                }
	            }

	            if (replacementUrl != null) {
	                resolvedUrl = resolvedUrl.replace("{" + placeholder + "}", replacementUrl);
	                usedExtensions.add(matchedExt);
	            }
	        }

	        // Remove query parameters that still contain unresolved placeholders
	        if (resolvedUrl.contains("{")) {
	            // Split URL into base and query string
	            String[] urlParts = resolvedUrl.split("\\?", 2);
	            if (urlParts.length == 2) {
	                String baseUrl = urlParts[0];
	                String queryString = urlParts[1];
	                
	                // Filter out parameters with placeholders
	                String[] params = queryString.split("&");
	                List<String> validParams = new ArrayList<>();
	                
	                for (String param : params) {
	                    if (!param.contains("{")) {
	                        validParams.add(param);
	                    }
	                }
	                
	                // Reconstruct URL
	                if (!validParams.isEmpty()) {
	                    resolvedUrl = baseUrl + "?" + String.join("&", validParams);
	                } else {
	                    resolvedUrl = baseUrl;
	                }
	            }
	        }

	        // Always add the URL, even if some placeholders were removed
	        if (!usedExtensions.isEmpty()) {
	            String enrichedLabel = "Send " + String.join(", ", usedExtensions) + " file(s) to " + toolLabel;
	            resolvedUrls.put(enrichedLabel, resolvedUrl);
	        } else {
	            // If no files were matched but URL is valid, add with generic label
	            resolvedUrls.put(toolLabel, resolvedUrl);
	        }
	    }
	    return resolvedUrls;
	}
}