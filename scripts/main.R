library(tidyverse)
library(ggplot2)
library(ape)
library(aplot)
library(geoR)
library(ggtree)
library(brotli)
library(brms)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

forms <- read.csv('../kinbank/cldf/forms.csv')
languages <- read.csv('../kinbank/cldf/languages.csv')

not_changing <- c('BornSameDay', 'HusbandsSameGroup', 'WifesSameGroup')

parameters <- read.csv('../kinbank/cldf/parameters.csv') %>% 
  rowwise() %>%
  mutate(ID_new = ifelse(grepl('exchange', ID) & !ID %in% not_changing, gsub('exchange', 'X', ID), ID),
         ID_new = ifelse(grepl('_not', ID_new) & !ID_new %in% not_changing, gsub('_not', '!', ID_new), ID_new),
         ID_new = ifelse(grepl('_or', ID_new) & !ID_new %in% not_changing, gsub('_or', '^', ID_new), ID_new),
         ID_new1 = ifelse(grepl('^m', ID_new) & !ID_new %in% not_changing, gsub('^m', 'mE', ID_new), ID_new),
         ID_new2 = ifelse(grepl('^f', ID_new1) & !ID_new1 %in% not_changing, gsub('^f', 'fE', ID_new1), ID_new1),
         ID_new3 = ifelse(grepl('S', ID_new2) & !ID_new2 %in% not_changing, gsub('S', 'mC', ID_new2), ID_new2),
         ID_new4 = ifelse(grepl('D', ID_new3) & !ID_new3 %in% not_changing, gsub('D', 'fC', ID_new3), ID_new3),
         ID_new5 = ifelse(grepl('B', ID_new4) & !ID_new4 %in% not_changing, gsub('B', 'mS', ID_new4), ID_new4),
         ID_new6 = ifelse(grepl('Z', ID_new5) & !ID_new5 %in% not_changing, gsub('Z', 'fS', ID_new5), ID_new5),
         ID_new7 = ifelse(grepl('F', ID_new6) & !ID_new6 %in% not_changing, gsub('F', 'mP', ID_new6), ID_new6),
         ID_new8 = ifelse(grepl('M', ID_new7) & !ID_new7 %in% not_changing, gsub('M', 'fP', ID_new7), ID_new7),
         ID_new9 = ifelse(grepl('H', ID_new8) & !ID_new8 %in% not_changing, gsub('H', 'mO', ID_new8), ID_new8),
         ID_new10 = ifelse(grepl('W', ID_new9) & !ID_new9 %in% not_changing, gsub('W', 'fO', ID_new9), ID_new9),
         revised_ID = ID_new10
  ) %>%
  select(-starts_with('ID_new'))

df <- forms %>% 
  left_join(languages %>% select(-Name), by=c('Language_ID' = 'ID')) %>% 
  left_join(parameters %>% mutate(ID = case_when(ID == 'mTW' ~ 'MT',
                                                 ID == 'fTW' ~ 'FT',
                                                 TRUE ~ ID)), by=c('Parameter_ID' = 'ID')) %>% 
  filter(Glottocode != 'aust1271',
         !(Source == 'BarthW2017' & Glottocode == 'stan1293'),
         !(Source == 'el-khaissi_2020' & Glottocode == 'stan1323'),
         !(Source == 'Evans2016' & Glottocode == 'russ1263'),
         !(Source == 'BarthW2017' & Glottocode == 'stan1295')) %>% 
  select(-Source) %>%
  distinct()

extensions <- df %>% 
  filter(Form != 'TRUE') %>%
  group_by(Glottolog_Name, Glottocode, Family, Latitude, Longitude, Form) %>% 
  distinct() %>%
  arrange(revised_ID) %>% 
  summarize(extension = list(revised_ID)) %>% 
  rowwise() %>%
  mutate(extension = paste(extension, collapse=','),
         #code = paste(as.character(memCompress(extension)), collapse = '')
  ) %>% ungroup()

extension_list_by_language <- extensions %>% 
  group_by(Glottolog_Name, Glottocode, Family, Latitude, Longitude) %>% 
  summarize(extension_list = list(extension)) %>% 
  rowwise() %>% 
  mutate(extension_list = paste(extension_list, collapse='/ '),
         compressed_code = paste(as.character(memCompress(extension_list)), collapse = ''),
         num_kinship_terms = str_count(pattern='/', string=extension_list) + 1,
         complexity = object.size(memCompress(extension_list))
  ) %>% 
  ungroup()

write.csv(extension_list_by_language, '../output_table/complexity_new.csv')

# plot

tree_orig <- read.nexus('../data/EDGE6635-merged-relabelled.tree')

#stolen code
to_keep <- tree_orig$tip.label %>% as.data.frame() %>% rename(tip.label = ".") %>% 
  separate(col = tip.label , into = c("Language_ID", "Name_EDGE"), remove = F, sep = 8) %>% 
  inner_join(languages, by=c('Language_ID'='Glottocode')) 

tree <- keep.tip(tree_orig, to_keep$tip.label) 
tree$tip.label <- to_keep$Language_ID %>% unique()


#tree$tip.label <- to_keep$Language_ID

p <- ggtree(tree)

p$data$label[p$data$isTip] <- to_keep$Language_ID %>% unique()


p2 <- ggplot(extension_list_by_language %>% filter(Glottocode %in% tree$tip.label), aes(x=complexity, y=Glottocode, color = Family)) +
  geom_point(alpha=0.6) +
  geom_text(aes(label=Glottolog_Name), size=1, nudge_x =0, color='black') +
  theme_classic(5) +
  theme(legend.position = 'none') +
  ylab('') +
  theme(axis.text.y=element_blank(),
        axis.ticks.y = element_blank())

pdf('../imgs/kinship_complexity_kinbank.pdf', width=10, height=35)
print(p2 %>% insert_left(p))
dev.off()


# analysis 
pcs <- read.csv('../data/compiled_table_20231213_pca.csv') %>% select(Glottocode, PC1) %>% distinct()

data_pc <- extension_list_by_language %>% left_join(pcs, by='Glottocode') %>%
  filter(!is.na(Glottocode),
         !is.na(Longitude),
         !is.na(Latitude),
         !is.na(PC1),
         Glottocode %in% tree$tip.label) %>%
  group_by(Glottocode) %>% slice_min(order_by=PC1, n=1) %>%
  mutate(Glottocode2 = Glottocode) %>%
  ungroup()

write.csv(data_pc, '../output_table/main.csv', row.names = FALSE)

# visualize
fit1 <- lm(complexity ~ PC1, data=data_pc)

pdf('../imgs/main.pdf', width=10, height=6)
ggplot(data_pc, aes(x=PC1, y=complexity)) +
  stat_smooth(method = "lm", col = "orange") +
  geom_point(alpha=0.6) +
  xlab('Sociopolitical Complexity Scale (SCI)') +
  ylab('Kinship System Complexity') +
  scale_x_reverse() +
  geom_label(aes(x = 4, y = 1200), hjust = 0, 
             label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5),
                           "\nIntercept =",signif(fit1$coef[[1]],5 ),
                           " \nSlope =",signif(fit1$coef[[2]], 5),
                           " \nP =",signif(summary(fit1)$coef[2,4], 5))) +
  theme_classic(15)
dev.off()

pdf('../imgs/main_by_family.pdf', width=15, height=10)
ggplot(data_pc, aes(x=PC1, y=complexity)) +
  stat_smooth(method = "lm", col = "orange", se=F) +
  facet_wrap(~Family) +
  geom_point(alpha=0.6) +
  xlab('Sociopolitical Complexity Scale (SCI)') +
  ylab('Kinship System Complexity') +
  scale_x_reverse() +
  # geom_label(aes(x = 4, y = 1200), hjust = 0, 
  #            label = paste("Adj R2 = ",signif(summary(fit1)$adj.r.squared, 5),
  #                          "\nIntercept =",signif(fit1$coef[[1]],5 ),
  #                          " \nSlope =",signif(fit1$coef[[2]], 5),
  #                          " \nP =",signif(summary(fit1)$coef[2,4], 5))) +
  theme_classic(15)
dev.off()

# analysis
kappa = 1 
phi_1 = c(1, 1.25) # "Local" version: (sigma, phi) First value is not used

A <- vcv.phylo(tree, corr=TRUE)

spatial_covar_mat_local = varcov.spatial(data_pc[,c("Longitude", "Latitude")], cov.pars = phi_1, kappa = kappa)$varcov
dimnames(spatial_covar_mat_local) = list(data_pc$Glottocode2, data_pc$Glottocode2)
spatial_covar_mat_local <- spatial_covar_mat_local / max(spatial_covar_mat_local)

set.seed(123)
model <- brm(data=data_pc,
             data2=list(A=A, spatial_covar_mat_local=spatial_covar_mat_local),
             family = 'gaussian',
             formula = complexity ~ PC1 + (1 | gr(Glottocode2, cov=spatial_covar_mat_local)) + (1 | gr(Glottocode, cov=A)),
             control = list(adapt_delta = 0.95),
             iter=10000,
             cores=4
)

summary(model)

# Results
#  Family: gaussian 
#   Links: mu = identity; sigma = identity 
# Formula: complexity ~ PC1 + (1 | gr(Glottocode2, cov = spatial_covar_mat_local)) + (1 | gr(Glottocode, cov = A)) 
#    Data: data_pc (Number of observations: 440) 
#   Draws: 4 chains, each with iter = 10000; warmup = 5000; thin = 1;
#          total post-warmup draws = 20000

# Group-Level Effects: 
# ~Glottocode (Number of levels: 437) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)   143.15     29.66    89.22   205.88 1.00     3823     7873

# ~Glottocode2 (Number of levels: 437) 
#               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sd(Intercept)    62.66     28.42     5.07   111.41 1.00      876     2356

# Population-Level Effects: 
#           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept   420.45     65.31   289.66   549.44 1.00     8243    10973
# PC1         -11.02      4.19   -19.24    -2.84 1.00    15207    13909

# Family Specific Parameters: 
#       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# sigma   138.80     10.72   116.71   157.81 1.00     1152     3104

# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).







# plot maps
world <- ne_countries(scale = 'medium', returnclass = 'sf')

cp <- tribble(
  ~Family, ~hex,
  "Afro-Asiatic", "#FFA500",
  "Anim",	"#00BFFF",
  "Atlantic-Congo",	"#D2691E",
  "Austroasiatic",	"#7FFF00",
  "Austronesian",	"#32CD32",
  "East Bird's Head",	"#7FFFD4",
  "Indo-European", 	"#DC143C",
  "Japonic",	"#FF69B4",
  "Kxa",	"#D2B48C",
  "Mayan",	"#DA70D6",
  "Nilotic",	"#F4A460",
  "North Halmahera",	"#00FFFF",
  "Nuclear Trans New Guinea",	"#0000FF",
  "Sino-Tibetan",	"#FF0000",
  "Surmic",	"#FFD700",
  "Tai-Kadai",	"#00B000",
  "Ticuna-Yuri",	"#9400D3",
  "Timor-Alor-Pantar",	"#1E90FF",
  "Turkic",	"#FFC0CB",
  "Algic",'#FC0FC0',
  "Arawakan", '#00AB68',
  "Baibai-Fas", '#F987C5',
  "Bilua", '#C921A4',
  "Border", "#E75480",
  "Bosavi", '#DA1884',
  "Cariban", '#C0DCE6',
  "Dravidian", '#41Fa9A',
  "Goilalan", '#Fa78DA',
  'Kolopom', '#FF9A8A',
  "Lower Sepik-Ramu", '#C291A4',
  "Mairasic", '#FE828C',
  'Maybrat-Karon', '#FC749C',
  'Nakh-Daghestanian', '#EEF44F',
  "Nimboranic", '#74b3ce',
  "Nuclear Torricelli", "#0abab5",
  "Pama-Nyungan", '#710c04',
  "Savosavo", '#03ac13',
  "Sentanic", '#32612d',
  "Sepik", "#74b72e",
  "Sko", '#5dbb63',
  "South Bougainville", '#00563b',
  "Taiap", "#b2ac88",
  'Tupian', '#005452',
  'Uralic', '#00072d',
  'Uto-Aztecan', '#4c0507',
  'Yawa-Saweru', '#5c004b',
  "Yele", '#500404',
  "", "#555555"
  )

sites <- st_as_sf(data_pc, coords = c("Longitude", "Latitude"), crs = 4326, agr = 'constant')  %>% left_join(cp, by='Family')

pdf('../imgs/map.pdf', width=10, height=6)
ggplot(data = world) +
  theme_bw() +
  geom_sf(linewidth=0.1) +
  geom_sf(data = sites, size = 2, shape = 23, fill = sites$hex, alpha=0.5, linewidth=0.1) +
  theme(legend.position = 'none')
dev.off()


# illustration of kinship complexity (toy examples)
ex <- tibble(
  scenarios = c(1,1,2,2,3,3),
  language = c('a','b','c','d','e','f'),
  extension_list = c(
    'fEfP',
    'fEfP/fEmP',
    'fEfP/fEfP/mEmP/mEmP',
    'fEfP/mEfP/mEmP/fEmP',
    'fEfSmC,fEmSmC,fEfSfC,fEmSfC',
    'mEmSfC,mEfPfS,fEfOmP,fEmPmS'
  )
) %>% 
  rowwise() %>%
  mutate(code=paste(as.character(memCompress(extension_list)), collapse = ''),
         length=str_count(code)) 


