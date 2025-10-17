import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

terrain_values = []
water_values_dict = {}

with open("./build/sampled_values.txt", "r") as file:
    lines = file.readlines()
    terrain_index = lines.index("Terrain sampled values:\n") + 1
    terrain_values = list(map(float, lines[terrain_index].strip().split()))

    water_lines = [line for line in lines if line.startswith("Water sampled values (timestamp")]
    for line in water_lines:
        timestamp = int(line.split("timestamp ")[1].split("):")[0])
        water_values = list(map(float, lines[lines.index(line) + 1].strip().split()))
        water_values_dict[timestamp] = water_values

plt.figure(figsize=(6, 4))
plt.plot(terrain_values, label="Terrain", color="green", alpha=0.7)

timestamps = list(water_values_dict.keys())
norm = plt.Normalize(min(timestamps), max(timestamps))


# colormap = cm.Blues
colormap = cm.winter
# colormap = cm.autumn
# colormap = cm.RdBu

for timestamp, water_values in water_values_dict.items():
    if timestamp % 5 == 0:
        filtered_water_values = [
            water if water >= terrain else None
            for water, terrain in zip(water_values, terrain_values)
        ]
        color = colormap(norm(timestamp))
        plt.plot(filtered_water_values, label=f"t={timestamp}", color=color, alpha=0.7)



# plt.title("Wetting and Drying Process")
plt.title("Tsunami Attenuation")
plt.xlabel("Timestamp")
plt.ylabel("Normalised Height")
# plt.legend()
plt.grid(alpha=0.3)
plt.savefig("wet_dry_plot.png", dpi=300, bbox_inches="tight")
plt.show()