# ==============================================================================
# Trait data curation
# ==============================================================================
# - Loads raw trait data and simplifies into binary (and 3-state, separating parasites)
#   categories.
# - Removes species with unknown/polymorphic/variable sociality/nesting behavior.
# - Saves curated trait table and a tip-matched, pruned tree to curated_data folder
# - Plots tip pies for each trait encoding (binary and three-state).
# ==============================================================================
# Setup
#-------------------------------------------------------------------------------
rm(list = ls())

data_dir <- "/home/lenarh/data/bee_nesting_sociality"

results_dir <- file.path(wd, "results", "supplemental")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

library(ape)
library(phytools)

#-------------------------------------------------------------------------------
# Load trait dataset and tree
# Prune trait dataset to just species in the tree and columns we will use
#-------------------------------------------------------------------------------
traits <- read.csv(
  file.path(wd, "original_data", "bee_traits.csv")
)

tree <- read.tree(
  file.path(wd, "original_data", "bee_tree_4586tips.tre")
)
Ntip(tree) # 4586
nrow(traits) # 4743 - contains some species not in current version of tree (from preprint version)

# Remove whitespace from species names first
tree$tip.label <- trimws(tree$tip.label)
traits$tips <- trimws(as.character(traits$tips))

# Find overlap
common_tips <- intersect(tree$tip.label, traits$tips) # 4586

# Trim trait dataset to only species matching the tree
traits <- traits[traits$tips %in% common_tips, ]

# Select and rename columns
traits <- traits[, c(
  "family",
  "tribe",
  "tips",
  "broad_sociality",
  "broad_nesting",
  "parasite_host_nesting"
)]

colnames(traits) <- c(
  "family",
  "tribe",
  "tips",
  "sociality",
  "nest",
  "parasite_nesting"
)

# Remove whitespace from remaining columns before manipulating
traits$sociality <- trimws(as.character(traits$sociality))
traits$nest <- trimws(as.character(traits$nest))
traits$parasite_nesting <- trimws(as.character(traits$parasite_nesting))

# Check data
table(traits$sociality)
# kleptoparasite     polymorphic     social     social parasite        solitary
#           461              36         972                 27            2983
# unknown
#     107

table(traits$nest)
# above-ground          ground  kleptoparasite  social parasite         unknown
#       1396              2570             461              27             122
# variable
# 10

head(traits)

# Now, we can recode sociality and nesting into binary categories.

#-------------------------------------------------------------------------------
# Sociality recode: Three-state (solitary, social, parasite)
#-------------------------------------------------------------------------------
# Initialize column
traits$sociality_three_states <- NA

# Create objects that describe what counts as solitary, social, parasite, or drop
solitary <- c("solitary")
social <- c("social")
drop <- c("unknown", "polymorphic")
parasite <- c("social parasite", "kleptoparasite")

# Assign "solitary" to all rows of sociality_binary where
# the sociality column matches any of the values in the solitary object
# and do the same for social, parasite, and drop
traits$sociality_three_states[which(
  traits$sociality %in% solitary
)] <- "solitary"
traits$sociality_three_states[which(traits$sociality %in% social)] <- "social"
traits$sociality_three_states[which(
  traits$sociality %in% parasite
)] <- "parasite"
traits$sociality_three_states[which(traits$sociality %in% drop)] <- "drop"

# Check counts
table(traits$sociality_three_states)
# drop parasite   social solitary
# 143      488      972     2983

#-------------------------------------------------------------------------------
# Sociality recode: Binary (solitary vs social)
#-------------------------------------------------------------------------------
traits$sociality_binary <- NA

solitary <- c("solitary", "kleptoparasite")
social <- c("social parasite", "social")
drop <- c("unknown", "polymorphic")


traits$sociality_binary[which(traits$sociality %in% solitary)] <- "solitary"
traits$sociality_binary[which(traits$sociality %in% social)] <- "social"
traits$sociality_binary[which(traits$sociality %in% drop)] <- "drop"

table(traits$sociality_binary)
# drop   social solitary
# 143      999     3444

#-------------------------------------------------------------------------------
# Remove “drop” rows for sociality (unknowns/polymorphic)
#-------------------------------------------------------------------------------
nrow(traits) # 4586
traits <- subset(traits, traits$sociality_binary != "drop")
nrow(traits) # 4443
traits <- subset(traits, traits$sociality_three_states != "drop")
nrow(traits) # 4443

4586 - 4443 # 143, same as the number of "drop"'s due to unknown/polymorphic sociality

#-------------------------------------------------------------------------------
# Nesting recode: Three-state (ground, aboveground, parasite)
#-------------------------------------------------------------------------------
traits$nest_three_states <- NA

ground <- c("ground")
aboveground <- c("above-ground")
drop <- c("unknown", "variable")
parasite <- c("social parasite", "kleptoparasite")

traits$nest_three_states[which(traits$nest %in% ground)] <- "ground"
traits$nest_three_states[which(traits$nest %in% aboveground)] <- "aboveground"
traits$nest_three_states[which(traits$nest %in% drop)] <- "drop"
traits$nest_three_states[which(traits$nest %in% parasite)] <- "parasite"

table(traits$nest_three_states)
# aboveground        drop      ground    parasite
#       1311         132        2512         488

#-------------------------------------------------------------------------------
# Nesting recode: Binary (ground vs aboveground)
#-------------------------------------------------------------------------------
# For parasites, replace their nesting state with that of their host
sum(
  traits$nest %in%
    c("kleptoparasite", "social parasite") & # Check how many
    traits$parasite_nesting %in% c("unknown", "variable"),
  na.rm = TRUE
) # 97

traits$nest[which(traits$nest == "kleptoparasite")] <- traits$parasite_nesting[
  traits$nest == "kleptoparasite"
]
traits$nest[which(traits$nest == "social parasite")] <- traits$parasite_nesting[
  traits$nest == "social parasite"
]

traits$nest_binary <- NA

ground <- c("ground")
aboveground <- c("above-ground")
drop <- c("unknown", "variable")

traits$nest_binary[which(traits$nest %in% ground)] <- "ground"
traits$nest_binary[which(traits$nest %in% aboveground)] <- "aboveground"
traits$nest_binary[which(traits$nest %in% drop)] <- "drop"

table(traits$nest_binary)
# aboveground        drop      ground
#       1330         229        2884

#-------------------------------------------------------------------------------
# Remove “drop” rows for nesting (unknown/variable)
#-------------------------------------------------------------------------------
nrow(traits) # 4443
traits <- subset(traits, traits$nest_binary != "drop")
nrow(traits) # 4214
traits <- subset(traits, traits$nest_three_states != "drop")
nrow(traits) # 4214

4443 - 4214 # 229, , same as the number of "drop"'s due to unknown/variable nesting

# So, the total number of species removed due to unknown/polymorphic/variable behavior is:
4586 - 4214 # 372

# Remove columns that are no longer needed (sociality, nest, parasite_nesting)
traits <- traits[, -c(4:6)]
colnames(traits)

#-------------------------------------------------------------------------------
# Prune tree to match trait dataset and save cleaned trait dataset and tree
#-------------------------------------------------------------------------------
Ntip(tree) # 4586
common_tips <- intersect(tree$tip.label, traits$tips)
tree <- keep.tip(tree, common_tips)
Ntip(tree) # 4214

write.tree(
  tree,
  file.path(wd, "curated_data", "bee_tree_pruned2.tre")
)

write.csv(
  traits,
  file.path(wd, "curated_data", "bee_traits_clean2.csv"),
  row.names = FALSE
)

#-------------------------------------------------------------------------------
# Visualize trait distributions on the pruned tree (tip pies)
#-------------------------------------------------------------------------------
# Align trait table to pruned tree tip order
colnames(traits)

group_traits <- traits[traits[, "tips"] %in% tree$tip.label, ]

# Order trait table to match tree tip order
group_traits <- group_traits[
  order(match(group_traits[, "tips"], tree$tip.label)),
]

#-------------------------------------------------------------------------------
# Sociality: three states
#-------------------------------------------------------------------------------
mode <- group_traits[, "sociality_three_states"]
names(mode) <- group_traits[, "tips"]

colors_states <- c("midnightblue", "goldenrod", "darkred")

pdf(
  file.path(results_dir, "sociality_three_states.pdf"),
  width = 4,
  height = 45
)

plot(tree, show.tip.label = TRUE, edge.width = 0.2, adj = 1, cex = 0.05)

par(fg = "transparent")

tiplabels(
  pie = to.matrix(mode, sort(unique(mode))),
  piecol = colors_states,
  cex = 0.1,
  lwd = 0.2,
  frame = "n"
)

par(fg = "black")

legend(
  "topleft",
  legend = sort(unique(mode)),
  pt.bg = colors_states,
  pch = 21,
  cex = 0.8
)

axisPhylo()
dev.off()

#-------------------------------------------------------------------------------
# Sociality: binary
#-------------------------------------------------------------------------------
mode <- group_traits[, "sociality_binary"]
names(mode) <- group_traits[, "tips"]

colors_states <- c("midnightblue", "goldenrod")

pdf(
  file.path(results_dir, "sociality_binary.pdf"),
  width = 4,
  height = 45
)

plot(tree, show.tip.label = TRUE, edge.width = 0.2, adj = 1, cex = 0.05)

par(fg = "transparent")

tiplabels(
  pie = to.matrix(mode, sort(unique(mode))),
  piecol = colors_states,
  cex = 0.1,
  lwd = 0.2,
  frame = "n"
)

par(fg = "black")

legend(
  "topleft",
  legend = sort(unique(mode)),
  pt.bg = colors_states,
  pch = 21,
  cex = 0.8
)

axisPhylo()
dev.off()

#-------------------------------------------------------------------------------
# Nesting: binary
#-------------------------------------------------------------------------------
mode <- group_traits[, "nest_binary"]
names(mode) <- group_traits[, "tips"]

colors_states <- c("midnightblue", "goldenrod")

pdf(
  file.path(results_dir, "nest_binary.pdf"),
  width = 4,
  height = 45
)

plot(tree, show.tip.label = TRUE, edge.width = 0.2, adj = 1, cex = 0.05)

par(fg = "transparent")

tiplabels(
  pie = to.matrix(mode, sort(unique(mode))),
  piecol = colors_states,
  cex = 0.1,
  lwd = 0.2,
  frame = "n"
)

par(fg = "black")

legend(
  "topleft",
  legend = sort(unique(mode)),
  pt.bg = colors_states,
  pch = 21,
  cex = 0.8
)

axisPhylo()
dev.off()

#-------------------------------------------------------------------------------
# Nesting: three states
#-------------------------------------------------------------------------------
mode <- group_traits[, "nest_three_states"]
names(mode) <- group_traits[, "tips"]

colors_states <- c("midnightblue", "goldenrod", "darkred")

pdf(
  file.path(results_dir, "nest_three_states.pdf"),
  width = 4,
  height = 45
)

plot(tree, show.tip.label = TRUE, edge.width = 0.2, adj = 1, cex = 0.05)

par(fg = "transparent")

tiplabels(
  pie = to.matrix(mode, sort(unique(mode))),
  piecol = colors_states,
  cex = 0.1,
  lwd = 0.2,
  frame = "n"
)

par(fg = "black")

legend(
  "topleft",
  legend = sort(unique(mode)),
  pt.bg = colors_states,
  pch = 21,
  cex = 0.8
)

axisPhylo()
dev.off()

#-------------------------------------------------------------------------------
# Done!
#-------------------------------------------------------------------------------
