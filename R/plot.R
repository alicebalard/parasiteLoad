# To keep all plot with the same theme (for Alice s TAC):
theme_alice <- theme_bw() +
  theme(plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 20),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        text = element_text(size = 20))

simdata

# Plot:
ggplot(Phuong_data, aes(x = HI, y = Oocysts_per_g)) +
  geom_point(size = 4) +
  geom_smooth() +
  theme_alice +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
