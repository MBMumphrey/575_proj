require(tidyverse)

human_metadata <- read_tsv(file="code/human.metadata")

normal <- filter(human_metadata,
                 diagnosis_s == "High Risk Normal" | diagnosis_s == "Normal") %>%
  sample_n(50)

adenoma <- filter(human_metadata,
                  diagnosis_s == "Adenoma" | diagnosis_s == "adv Adenoma") %>%
  sample_n(50)

cancer <- filter(human_metadata,
                 diagnosis_s == "Cancer") %>%
  sample_n(50)

subsample <- rbind(normal, adenoma, cancer)

human.files <- tibble(sample = subsample$Sample_Name_s,
                      fastq = subsample$Run_s) %>%
  mutate(fastq1 = paste(fastq, "_1.fastq", sep = "")) %>%
  mutate(fastq2 = paste(fastq, "_2.fastq", sep = "")) %>%
  select(sample, fastq1, fastq2)

write_tsv(human.files, "data/human.files", col_names = FALSE)
write_tsv(subsample, "data/human.metadata.subsample")
