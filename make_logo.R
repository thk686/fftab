library(hexSticker)
library(RColorBrewer)
library(ggplot2)

n <- 5
m <- 200

t <- seq(-pi, pi, len = m)

fb <- list()
for (i in 1:n) {
  fb[[length(fb) + 1]] <- sin(i * t) / i - (i / 2)
  fb[[length(fb) + 1]] <- cos(i * t) / i - (i / 2)
}

df <- tibble::tibble(g = as.factor(rep(1:(2 * n), each = m)),
                     x = rep(t, (2 * n)), y = unlist(fb))

p <- ggplot(df) +
  geom_line(aes(x = x, y = y, color = g), lwd = 1) +
  scale_color_brewer(palette = "Paired") +
  guides(color = "none") +
  theme_void() +
  theme_transparent()

sticker(p,
        package="fftab",
        p_y = 1.63,
        p_size = 23,
        p_color = "black",
        s_x=1,
        s_y=0.95,
        s_width=2,
        s_height=1,
        h_size = 1,
        h_fill = "lightgrey",
        h_color = "black",
        filename="man/figures/logo.png")
