# sizeclass: divides the position measure for each FG by the total number of times that 
# FG appears in any position within the same motif size class (the number of nodes a motif 
# contains).
# "sizeclass_NAzero": same as ’sizeclass’ but replaces all NA values with 0. If a FG
# does not occur in any motifs in a given size class, ’sizeclass’ normalisation will return
# NAs. ’sizeclass_NAzero’ avoids this by replacing NAs with zero.

sizeclass_normalization <- function(positions_i_agg_raw){
  
  normalization <- positions_i_agg_raw %>% mutate(
    poll2node = np2,
    poll3node = np4 + np6,
    poll4node = np8 + np11 + np12 + np14 + np16,
    poll5node = np18 + np21 + np22 + np24 + np25 + np28 + np29 + np31 + np34 + np35 + np38 + 
      np41 + np42 + np44 + np46,
    plant2node = np1,
    plant3node = np3 + np5,
    plant4node = np7 + np9 + np10 + np13 + np15,
    plant5node = np17 + np19 + np20 + np23 + np26 + np27 + np30 + np32 + np33 + np36 + np37 + 
      np39 + np40 + np43 + np45,
    np1 = np1 / plant2node,
    np2 = np2 / poll2node,
    np3 = np3 / plant3node,
    np4 = np4 / poll3node,
    np5 = np5 / plant3node,
    np6 = np6 / poll3node,
    np7 = np7 / plant4node,
    np8 = np8 / poll4node,
    np9 = np9 / plant4node,
    np10 = np10 / plant4node,
    np11 = np11 / poll4node,
    np12 = np12 / poll4node,
    np13 = np13 / plant4node,
    np14 = np14 / poll4node,
    np15 = np15 / plant4node,
    np16 = np16 / poll4node,
    np17 = np17 / plant5node,
    np18 = np18 / poll5node,
    np19 = np19 / plant5node,
    np20 = np20 / plant5node,
    np21 = np21 / poll5node,
    np22 = np22 / poll5node,
    np23 = np23 / plant5node,
    np24 = np24 / poll5node,
    np25 = np25 / poll5node,
    np26 = np26 / plant5node,
    np27 = np27 / plant5node,
    np28 = np28 / poll5node,
    np29 = np29 / poll5node,
    np30 = np30 / plant5node,
    np31 = np31 / poll5node,
    np32 = np32 / plant5node,
    np33 = np33 / plant5node,
    np34 = np34 / poll5node,
    np35 = np35 / poll5node,
    np36 = np36 / plant5node,
    np37 = np37 / plant5node,
    np38 = np38 / poll5node,
    np39 = np39 / plant5node,
    np40 = np40 / plant5node,
    np41 = np41 / poll5node,
    np42 = np42 / poll5node,
    np43 = np43 / plant5node,
    np44 = np44 / poll5node,
    np45 = np45 / plant5node,
    np46 = np46 / poll5node
  ) %>% select(-poll2node, -poll3node, -poll4node, -poll5node,
               -plant2node, -plant3node, -plant4node, -plant5node) %>% ungroup()
  
  normalization[is.na(normalization)] <- 0
  
  return(normalization)
  
} 