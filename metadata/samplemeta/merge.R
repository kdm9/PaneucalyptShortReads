#install.packages(c("tidyverse", "stringdist"))
library(tidyverse)
library(stringdist)
library(lubridate)


meta = read_csv("../sample2runlib.csv") %>%
    filter(Include=="Y") %>%
    select(SampleID=sample) %>%
    unique() %>%
    as.tibble()


dp15 = read_csv("dp15.csv") %>%
    select(SampleID=ID, Species=Species, Location, Date, Latitude, Longitude, Elevation) %>%
    mutate(SampleName=SampleID, Species=sprintf("Eucalyptus %s", Species), Date=as.POSIXct(Date)) %>%
	mutate(Datum="WGS84") %>%
    inner_join(meta, by="SampleID")

lbmel = read_csv("LBEmel_BGI.csv") %>%
    filter(HaveBGISeq) %>%
    mutate(Altitude=as.numeric(Altitude), Species="Eucalyptus melliodora") %>%
    select(SampleID=AnonName, SampleName=AccessionName, Species, Latitude, Longitude, Elevation=Altitude, Population=PopulationName) %>% 
	mutate(Datum="WGS84") %>%
    inner_join(meta, by="SampleID")

cca = read_tsv("CCACurrentWithMeta.csv") %>%
    select(SampleID=accession_id, Species=CurrentName, Location, Latitude, Longitude) %>% 
	mutate(Datum="WGS84") %>%
    unique() %>%
    inner_join(meta, by="SampleID")


pan = read_csv("pangenome.csv") %>%
	mutate(Datum="WGS84")

hmb = read_csv("HMB2019Collections_All.csv") %>%
	mutate(Date=dmy_hms(Date)) %>%
	select(SampleID=Ident, Species, Datum=Coord_system, Latitude, Longitude, Date, Elevation=Altitude_m, Location) %>%
	inner_join(meta, by="SampleID")


filled_meta = bind_rows(dp15, lbmel, cca, pan, hmb) %>%
    left_join(meta, ., by="SampleID")

dups = table(filled_meta$SampleID)
dups[dups>1]

write_csv(filled_meta, "../sample-metadata.csv")

