library(tidyverse)

samp = read_csv("sample-metadata.csv")
runs = read_csv("sample2runlib.csv")

proj = runs %>%
    select(Sample=sample, Project=project) %>%
    unique()

if (!all(table(proj$Sample) == 1))
    stop("Not all samples have only one project")

projsum = proj %>% 
    left_join(samp, by=c("Sample"="SampleID")) %>%
    group_by(Project) %>%
    summarise(n_indiv=n(), n_spp=length(unique(Species)))
write_csv(projsum, "summaries/project.csv")

sppsum = samp %>%
    filter(!is.na(Species)) %>%
    group_by(Species) %>%
    summarise(n_indiv=n(), n_loc=length(unique(paste(Latitude, Longitude)))) %>%
    ungroup() %>%
    arrange(-n_indiv)
write_csv(sppsum, "summaries/species.csv")
