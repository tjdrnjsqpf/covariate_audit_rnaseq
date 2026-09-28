#!/usr/bin/env Rscript
## Print at once all content-of-change (46) numbers that enter the v4 manuscript (based on C_content_by_class.tsv). Usage: Rscript scripts/06_revision/51_v4_numbers.R [suffix]
suppressMessages(library(data.table)); suf <- commandArgs(trailingOnly=TRUE)[1]; if (is.na(suf)) suf <- ""
S <- fread(sprintf("results/revision/C_content_by_class%s.tsv", suf)); lv <- c("balanced","mildly_imbalanced","partially_confounded")
g <- function(o, k, col) S[organism == o & design_sex == k][[col]]
pc <- function(x) sprintf("%.0f%%", 100*as.numeric(x)); m3 <- function(x) sub(" \\[.*", "", x)
cat("n by class:", paste(S[, sprintf("%s %s=%d", organism, design_sex, n)], collapse="; "), "\n")
cat("B:", paste(S[, sprintf("%s %s=%g", organism, design_sex, B)], collapse="; "), "\n\n")
k <- "partially_confounded"
cat("== Results §5 ¶3 (confounded) ==\n")
cat("below all within perms (k=0):", pc(g("human",k,"frac_p_within_le05")), pc(g("mouse",k,"frac_p_within_le05")), " | below all across perms:", pc(g("human",k,"frac_p_across_le05")), pc(g("mouse",k,"frac_p_across_le05")), "\n")
cat("lost overlap real-vs-within:", g("human",k,"lost_overlap_real_vs_within"), "|", g("mouse",k,"lost_overlap_real_vs_within"), "\n")
cat("lost overlap null:", g("human",k,"lost_overlap_null"), "|", g("mouse",k,"lost_overlap_null"), "\n")
cat("ratio:", g("human",k,"lost_overlap_ratio"), "|", g("mouse",k,"lost_overlap_ratio"), "\n")
cat("xy_lost real/within:", g("human",k,"xy_lost_M1"), g("human",k,"xy_lost_within"), "|", g("mouse",k,"xy_lost_M1"), g("mouse",k,"xy_lost_within"), "\n")
cat("d_within:", g("human",k,"d_within"), "|", g("mouse",k,"d_within"), "   d_within_noxy:", g("human",k,"d_within_noxy"), "|", g("mouse",k,"d_within_noxy"), "\n")
cat("lfc shift real/within:", g("human",k,"lfc_shift_M1"), g("human",k,"lfc_shift_within"), "|", g("mouse",k,"lfc_shift_M1"), g("mouse",k,"lfc_shift_within"), "\n")
cat("xy share of shift real/within:", g("human",k,"xy_share_shift_M1"), g("human",k,"xy_share_shift_within"), "|", g("mouse",k,"xy_share_shift_M1"), g("mouse",k,"xy_share_shift_within"), "\n")
cat("tsp real/within:", g("human",k,"tsp_M1"), g("human",k,"tsp_within"), "|", g("mouse",k,"tsp_M1"), g("mouse",k,"tsp_within"), "\n")
cat("\n== balanced (Fig 4A sentence) ==\n"); k2 <- "balanced"
cat("below all within perms (k=0):", pc(g("human",k2,"frac_p_within_le05")), pc(g("mouse",k2,"frac_p_within_le05")), " | below all across perms:", pc(g("human",k2,"frac_p_across_le05")), pc(g("mouse",k2,"frac_p_across_le05")), " (null floor 1/(B+1) =", sprintf("%.0f%%", 100/(g("human",k2,"B")+1)), ")\n")
cat("d_within balanced:", g("human",k2,"d_within"), "|", g("mouse",k2,"d_within"), "\n")
cat("\n== Results §7 (pathway, confounded) ==\n")
cat("PJ below all within perms:", pc(g("human",k,"frac_p_PJ_within_le05")), pc(g("mouse",k,"frac_p_PJ_within_le05")), " | p_across:", pc(g("human",k,"frac_p_PJ_across_le05")), pc(g("mouse",k,"frac_p_PJ_across_le05")), "\n")
cat("top1 changed real / within(mean over perms) / across:", pc(g("human",k,"top1_changed_M1")), pc(g("human",k,"top1_changed_within")), pc(g("human",k,"top1_changed_across")), "|", pc(g("mouse",k,"top1_changed_M1")), pc(g("mouse",k,"top1_changed_within")), pc(g("mouse",k,"top1_changed_across")), "\n")
