library(tidyverse)

samp = read_csv("sample-metadata.csv")

x = data.frame(table(samp$Species))
colnames(x) = c("Species", "N")
write_csv(x, "species.csv")

spp = read_csv("species.csv")

sets = samp %>%
    inner_join(spp, by="Species")

write_csv(sets, "samples-by-sppref.csv")


x = sets %>%
    group_by(SetName) %>%
    group_walk(function (tbl, setname) {
           tbl %>% pull(SampleID) %>% write_lines(sprintf("samplesets/SPP_%s.txt", setname))
    })
