# ============================================================================
# render_warmup_mechanismB_flowchart.R
# Offline static PNG of the Mechanism B (profileInitOnStall) on-stall rescue,
# drawn natively with ggplot2 so it needs no browser / Graphviz / network.
# Companion to presentation/_diagrams/profileInit_warmup_mechanismB.md.
#
# Usage:  source("script/viz/render_warmup_mechanismB_flowchart.R")
# ============================================================================

suppressPackageStartupMessages({library(ggplot2)})

OUT <- "presentation/_diagrams/profileInit_warmup_mechanismB.png"

## node table: id, x, y, half-width, half-height, type, label
nd <- function(id, x, y, w, h, type, label)
  data.frame(id, x, y, w, h, type, label, stringsAsFactors = FALSE)

nodes <- rbind(
  nd("A", 0,  9, 3.3, 0.7, "io",   "Forward candidate (var ~ covar . shape), add == TRUE"),
  nd("B", 0,  8, 3.0, 0.6, "proc", "Fit candidate with real estimator (focei / bobyqa)"),
  nd("C", 0,  7, 3.0, 0.6, "proc", "dOFV = fit$objf - candidate$objf"),
  nd("D", 0, 5.9, 3.2, 0.9, "dec", "profileInitOnStall & add &\nfinite(dOFV) & dOFV <= stallTol ?"),
  nd("Z", 4.4, 5.9, 2.3, 0.7, "term", "Keep candidate as-is\n(healthy: dOFV > stallTol)"),
  nd("E", 0, 4.6, 3.4, 0.9, "profl", ".profileCovInit(): 1-D FOCEi profile\n(BOBYQA outer, others fixed, BSV fixed-present,\nEBEs re-estimated each trial)"),
  nd("F", 0, 3.3, 3.2, 0.8, "dec", "Finite & strictly inside\n[lower+tol, upper-tol] ?"),
  nd("G", 4.4, 3.3, 2.3, 0.6, "term", "Return NA -> discard rescue"),
  nd("H", 0, 2.1, 3.2, 0.7, "proc", "Refit real estimator from profiled init -> x_p"),
  nd("I", 0, 1.2, 3.0, 0.55, "proc", "dOFV_p = fit$objf - x_p$objf"),
  nd("J", 0, 0.1, 3.0, 0.8, "dec", "x_p valid &\ndOFV_p > dOFV ?"),
  nd("K", 4.4, 0.1, 2.3, 0.6, "term", "Reject rescue\n(keep stalled candidate)"),
  nd("L", 0, -1.2, 3.2, 0.7, "accept", "Accept rescue: replace x, dOFV,\ndof, pchisqr with x_p"),
  nd("M", 0, -2.3, 3.4, 0.6, "io", "Record step-table row -> selection")
)

## edges: from-id, to-id, label, side ("v" straight down, or route via elbow)
ed <- function(from, to, label = "", type = "v")
  data.frame(from, to, label, type, stringsAsFactors = FALSE)

edges <- rbind(
  ed("A","B"), ed("B","C"), ed("C","D"),
  ed("D","Z","No"), ed("D","E","Yes (stalled)"),
  ed("E","F"),
  ed("F","G","No"), ed("F","H","Yes"),
  ed("H","I"), ed("I","J"),
  ed("J","K","No"), ed("J","L","Yes (better)"),
  ed("L","M")
)

pal <- c(io = "#E8EAED", proc = "#D6E4F0", dec = "#FCE5B6",
         profl = "#D8F0DC", accept = "#CDEBD6", term = "#F0D9D9")

pos <- function(id) nodes[match(id, nodes$id), ]

## build edge segments (elbow for horizontal targets)
seg <- do.call(rbind, lapply(seq_len(nrow(edges)), function(i) {
  e <- edges[i, ]; a <- pos(e$from); b <- pos(e$to)
  if (abs(a$x - b$x) < 1e-6) {           # vertical
    data.frame(x = a$x, y = a$y - a$h, xend = b$x, yend = b$y + b$h,
               lx = a$x + 0.12, ly = (a$y - a$h + b$y + b$h)/2, lab = e$label)
  } else {                                # elbow to the right
    data.frame(x = a$x + a$w, y = a$y, xend = b$x - b$w, yend = b$y,
               lx = (a$x + a$w + b$x - b$w)/2, ly = b$y + 0.22, lab = e$label)
  }
}))

p <- ggplot() +
  geom_segment(data = seg, aes(x, y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
               linewidth = 0.4, colour = "grey30") +
  geom_text(data = seg, aes(lx, ly, label = lab),
            size = 2.5, hjust = 0, colour = "grey20") +
  geom_tile(data = nodes, aes(x, y, width = 2*w, height = 2*h, fill = type),
            colour = "grey40", linewidth = 0.4) +
  geom_text(data = nodes, aes(x, y, label = label), size = 2.5, lineheight = 0.9) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "SCM warm-start - Mechanism B (profileInitOnStall): on-stall rescue") +
  coord_equal(clip = "off") +
  theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.margin = margin(10, 10, 10, 10))

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
ggsave(OUT, p, width = 9, height = 12, dpi = 200, bg = "white")
message("saved: ", OUT)
