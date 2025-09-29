library(httr)
library(jsonlite)
library(tidyverse)

########## Entity Schema #############
# Make GET request to the URL
response <- GET("https://data.rcsb.org/rest/v1/schema/entry")

# Extract content
content <- content(response, "text")

# Parse JSON content into an R list
parsed_content <- fromJSON(content)

# Initialize an empty data frame to store properties and their sub-properties
properties_df <- as.data.frame(matrix(data = NA, byrow = length(names(parsed_content$properties)), ncol = 2))
colnames(properties_df) = c("Category", "Properties")
head(properties_df)

# Loop through each property to fetch sub-properties
for (i in 1:length(names(parsed_content$properties))) {

  print(i)
  prop = names(parsed_content$properties)[i]

  sub_props = names(parsed_content$properties[[prop]]$items$properties)
  properties_df[i,1] = prop


  if(is.null(sub_props)){

    sub_props = names(parsed_content$properties[[prop]]$properties)

    properties_df[i,2] = paste(sub_props, collapse = ", ")



  }else{

    properties_df[i,2] = paste(sub_props, collapse = ", ")

  }

}

library(writexl)
# Write the results to a CSV file
write_xlsx(properties_df, "RCSB_properties_and_sub_properties_entry.xlsx")


########## Polymer Entity Schema #############
# Make GET request to the URL
response <- GET("https://data.rcsb.org/rest/v1/schema/polymer_entity")

# Extract content
content <- content(response, "text")

# Parse JSON content into an R list
parsed_content <- fromJSON(content)

# Initialize an empty data frame to store properties and their sub-properties
properties_df <- as.data.frame(matrix(data = NA, byrow = length(names(parsed_content$properties)), ncol = 2))
colnames(properties_df) = c("Category", "Properties")
head(properties_df)

# Loop through each property to fetch sub-properties
for (i in 1:length(names(parsed_content$properties))) {

  print(i)
  prop = names(parsed_content$properties)[i]

  sub_props = names(parsed_content$properties[[prop]]$items$properties)
  properties_df[i,1] = prop


  if(is.null(sub_props)){

    sub_props = names(parsed_content$properties[[prop]]$properties)

    properties_df[i,2] = paste(sub_props, collapse = ", ")



  }else{

    properties_df[i,2] = paste(sub_props, collapse = ", ")

  }

}

library(writexl)
# Write the results to a CSV file
write_xlsx(properties_df, "data/RCSB_properties_and_sub_properties_polymer_entity.xlsx")



########## Branched Entity Schema #############
# Make GET request to the URL
response <- GET("https://data.rcsb.org/rest/v1/schema/branched_entity")

# Extract content
content <- content(response, "text")

# Parse JSON content into an R list
parsed_content <- fromJSON(content)

# Initialize an empty data frame to store properties and their sub-properties
properties_df <- as.data.frame(matrix(data = NA, byrow = length(names(parsed_content$properties)), ncol = 2))
colnames(properties_df) = c("Category", "Properties")
head(properties_df)

# Loop through each property to fetch sub-properties
for (i in 1:length(names(parsed_content$properties))) {

  print(i)
  prop = names(parsed_content$properties)[i]

  sub_props = names(parsed_content$properties[[prop]]$items$properties)
  properties_df[i,1] = prop


  if(is.null(sub_props)){

    sub_props = names(parsed_content$properties[[prop]]$properties)

    properties_df[i,2] = paste(sub_props, collapse = ", ")



  }else{

    properties_df[i,2] = paste(sub_props, collapse = ", ")

  }

}

library(writexl)
# Write the results to a CSV file
write_xlsx(properties_df, "data/RCSB_properties_branched_entity.xlsx")


########## Non Polymer Entity Schema #############
# Make GET request to the URL
response <- GET("https://data.rcsb.org/rest/v1/schema/nonpolymer_entity")

# Extract content
content <- content(response, "text")

# Parse JSON content into an R list
parsed_content <- fromJSON(content)

# Initialize an empty data frame to store properties and their sub-properties
properties_df <- as.data.frame(matrix(data = NA, byrow = length(names(parsed_content$properties)), ncol = 2))
colnames(properties_df) = c("Category", "Properties")
head(properties_df)

# Loop through each property to fetch sub-properties
for (i in 1:length(names(parsed_content$properties))) {

  print(i)
  prop = names(parsed_content$properties)[i]

  sub_props = names(parsed_content$properties[[prop]]$items$properties)
  properties_df[i,1] = prop


  if(is.null(sub_props)){

    sub_props = names(parsed_content$properties[[prop]]$properties)

    properties_df[i,2] = paste(sub_props, collapse = ", ")



  }else{

    properties_df[i,2] = paste(sub_props, collapse = ", ")

  }

}

library(writexl)
# Write the results to a CSV file
write_xlsx(properties_df, "data/RCSB_properties_nonpolymer_entity.xlsx")


########## Polymer Instance Schema #############
# Make GET request to the URL
response <- GET("https://data.rcsb.org/rest/v1/schema/polymer_entity_instance")

# Extract content
content <- content(response, "text")

# Parse JSON content into an R list
parsed_content <- fromJSON(content)

# Initialize an empty data frame to store properties and their sub-properties
properties_df <- as.data.frame(matrix(data = NA, byrow = length(names(parsed_content$properties)), ncol = 2))
colnames(properties_df) = c("Category", "Properties")
head(properties_df)

# Loop through each property to fetch sub-properties
for (i in 1:length(names(parsed_content$properties))) {

  print(i)
  prop = names(parsed_content$properties)[i]

  sub_props = names(parsed_content$properties[[prop]]$items$properties)
  properties_df[i,1] = prop


  if(is.null(sub_props)){

    sub_props = names(parsed_content$properties[[prop]]$properties)

    properties_df[i,2] = paste(sub_props, collapse = ", ")



  }else{

    properties_df[i,2] = paste(sub_props, collapse = ", ")

  }

}

library(writexl)
# Write the results to a CSV file
write_xlsx(properties_df, "data/RCSB_properties_polymer_entity_instance.xlsx")


########## Branched Instance Schema #############
# Make GET request to the URL
response <- GET("https://data.rcsb.org/rest/v1/schema/branched_entity_instance")

# Extract content
content <- content(response, "text")

# Parse JSON content into an R list
parsed_content <- fromJSON(content)

# Initialize an empty data frame to store properties and their sub-properties
properties_df <- as.data.frame(matrix(data = NA, byrow = length(names(parsed_content$properties)), ncol = 2))
colnames(properties_df) = c("Category", "Properties")
head(properties_df)

# Loop through each property to fetch sub-properties
for (i in 1:length(names(parsed_content$properties))) {

  print(i)
  prop = names(parsed_content$properties)[i]

  sub_props = names(parsed_content$properties[[prop]]$items$properties)
  properties_df[i,1] = prop


  if(is.null(sub_props)){

    sub_props = names(parsed_content$properties[[prop]]$properties)

    properties_df[i,2] = paste(sub_props, collapse = ", ")



  }else{

    properties_df[i,2] = paste(sub_props, collapse = ", ")

  }

}

library(writexl)
# Write the results to a CSV file
write_xlsx(properties_df, "data/RCSB_properties_branched_entity_instance.xlsx")


########## Non-polymer Instance Schema #############
# Make GET request to the URL
response <- GET("https://data.rcsb.org/rest/v1/schema/nonpolymer_entity_instance")

# Extract content
content <- content(response, "text")

# Parse JSON content into an R list
parsed_content <- fromJSON(content)

# Initialize an empty data frame to store properties and their sub-properties
properties_df <- as.data.frame(matrix(data = NA, byrow = length(names(parsed_content$properties)), ncol = 2))
colnames(properties_df) = c("Category", "Properties")
head(properties_df)

# Loop through each property to fetch sub-properties
for (i in 1:length(names(parsed_content$properties))) {

  print(i)
  prop = names(parsed_content$properties)[i]

  sub_props = names(parsed_content$properties[[prop]]$items$properties)
  properties_df[i,1] = prop


  if(is.null(sub_props)){

    sub_props = names(parsed_content$properties[[prop]]$properties)

    properties_df[i,2] = paste(sub_props, collapse = ", ")



  }else{

    properties_df[i,2] = paste(sub_props, collapse = ", ")

  }

}

library(writexl)
# Write the results to a CSV file
write_xlsx(properties_df, "data/RCSB_properties_nonpolymer_entity_instance.xlsx")


########## Assembly Schema #############
# Make GET request to the URL
response <- GET("https://data.rcsb.org/rest/v1/schema/assembly")

# Extract content
content <- content(response, "text")

# Parse JSON content into an R list
parsed_content <- fromJSON(content)

# Initialize an empty data frame to store properties and their sub-properties
properties_df <- as.data.frame(matrix(data = NA, byrow = length(names(parsed_content$properties)), ncol = 2))
colnames(properties_df) = c("Category", "Properties")
head(properties_df)

# Loop through each property to fetch sub-properties
for (i in 1:length(names(parsed_content$properties))) {

  print(i)
  prop = names(parsed_content$properties)[i]

  sub_props = names(parsed_content$properties[[prop]]$items$properties)
  properties_df[i,1] = prop


  if(is.null(sub_props)){

    sub_props = names(parsed_content$properties[[prop]]$properties)

    properties_df[i,2] = paste(sub_props, collapse = ", ")



  }else{

    properties_df[i,2] = paste(sub_props, collapse = ", ")

  }

}

library(writexl)
# Write the results to a CSV file
write_xlsx(properties_df, "data/RCSB_properties_assembly.xlsx")


########## Chemical Component Schema #############
# Make GET request to the URL
response <- GET("https://data.rcsb.org/rest/v1/schema/chem_comp")

# Extract content
content <- content(response, "text")

# Parse JSON content into an R list
parsed_content <- fromJSON(content)

# Initialize an empty data frame to store properties and their sub-properties
properties_df <- as.data.frame(matrix(data = NA, byrow = length(names(parsed_content$properties)), ncol = 2))
colnames(properties_df) = c("Category", "Properties")
head(properties_df)

# Loop through each property to fetch sub-properties
for (i in 1:length(names(parsed_content$properties))) {

  print(i)
  prop = names(parsed_content$properties)[i]

  sub_props = names(parsed_content$properties[[prop]]$items$properties)
  properties_df[i,1] = prop


  if(is.null(sub_props)){

    sub_props = names(parsed_content$properties[[prop]]$properties)

    properties_df[i,2] = paste(sub_props, collapse = ", ")



  }else{

    properties_df[i,2] = paste(sub_props, collapse = ", ")

  }

}

library(writexl)
# Write the results to a CSV file
write_xlsx(properties_df, "data/RCSB_properties_chem_comp.xlsx")
