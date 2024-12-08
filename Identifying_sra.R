
#Identifying the right SRA entries

sra_search <- entrez_search(db = "sra", term = "Bombus terrestris", use_history = TRUE, retmax = 6000) 
length(sra_search$ids)  
sra_ids_table <- data.frame(SRR_IDs = sra_search$ids)
print(sra_ids_table)

#saving SRA IDS as a table 
write.table(sra_ids_table, file = "sra_ids.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


# Use BioProject database instead of SRA to search for Project IDs (PRJNA)
search_term <- "Bombus terrestris"
bioproject_search <- entrez_search(db = "bioproject", term = search_term, use_history = TRUE, retmax = 6000)

print(paste("Number of Project IDs found:", length(bioproject_search$ids)))

project_ids_table <- data.frame(Project_IDs = bioproject_search$ids)

print(project_ids_table)

# write the Project IDs to a txt file
write.table(project_ids_table, file = "C:/Users/User/Documents/Masterthesis/bioproject_ids.txt", row.names = FALSE)