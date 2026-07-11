# US Study Sites Map
# Plots a map of the US with Nebraska, Alabama, and Georgia highlighted
# and points marking Lincoln, NE; Huntsville, AL; and Fort Valley, GA

library(ggplot2)
library(maps)
library(paletteer)

# Get US state map data
us_states <- map_data("state")

# Exclude states west of ND/SD/NE/KS/OK/TX
western_states <- c("montana", "wyoming", "colorado", "new mexico",
                    "idaho", "utah", "arizona", "nevada", "california",
                    "oregon", "washington", "hawaii")
us_states <- us_states[!us_states$region %in% western_states, ]

# Define states to highlight
highlight_states <- c("nebraska", "alabama", "georgia")
# n_genotypes <- c(972, 340, 692)
# n_images <- c(6147, 1735, 5634)

# Create a column to indicate which states should be highlighted
us_states$highlight <- ifelse(us_states$region %in% highlight_states, "yes", "no")
# us_states <- us_states %>% 
#   mutate(n_genotypes = case_when(region=='nebraska' ~ 972, 
#                                  region=='georgia' ~ 692, 
#                                  region=='alabama' ~ 340, 
                                 # .default = NA))

# City coordinates
cities <- data.frame(
  city = c("Lincoln", "Huntsville", "Fort Valley"),
  state = c("Nebraska", "Alabama", "Georgia"),
  lat = c(40.861, 34.902, 33.561),
  lon = c(-96.597, -86.561, -84.313), 
  n_images = c(6147, 1735, 3695)
)
# pal_badlands <- paletteer_d('nationalparkcolors::Badlands', 5)
# Create the map
p <- ggplot() +
  # Draw all states
  geom_polygon(data = us_states,
               aes(x = long, y = lat, group = group, fill = highlight),
               color = "black", linewidth = 0.2) +
  # Set fill colors
  scale_fill_manual(values = c("yes" = paletteer_d("nationalparkcolors::Acadia")[3], "no" = "white"),
                    guide = "none") +
  # scale_fill_viridis(option = 'plasma', direction = -1, name = 'Genotypes') +
  # Add city points
  geom_point(data = cities,
             aes(x = lon, y = lat),
             color = 'black', size = 5) +
  # Add city labels
  # geom_text(data = cities,
  #           aes(x = lon, y = lat, label = paste(city, state, sep = ", ")),
  #           hjust = -0.1, vjust = 0.5, size = 3) +
  # Use map projection

  coord_map("albers", lat0 = 29.5, lat1 = 45.5) +
  # Clean theme
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

# Display the plot
print(p)

# Save the plot
ggsave("us_study_sites_map.png", p, width = 6.4, height = 5.2, dpi = 300, bg = "white")
