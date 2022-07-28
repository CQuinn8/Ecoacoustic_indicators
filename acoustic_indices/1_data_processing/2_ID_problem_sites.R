
library(dplyr)

df = read.csv('//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/data/Recordings-Grid view.csv')

problem_sites_df = df %>% 
  select(UniqueID, SiteID, Notes) %>%
  filter(!Notes == "")

write.csv(problem_sites_df, 
          '//shares.hpc.nau.edu/cirrus/projects/tropics/users/cquinn/s2l/paper1-AcousticIndices/data/identified_problem_sites.csv',
          row.names = FALSE)
