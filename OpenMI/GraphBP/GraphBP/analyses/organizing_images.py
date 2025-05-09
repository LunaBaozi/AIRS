import math
import os
import matplotlib.pyplot as plt

# Config:
images_dir = './your_dir_with_images'
result_grid_filename = './grid.jpg'
result_figsize_resolution = 40 # 1 = 100px

images_list = os.listdir(images_dir)
images_count = len(images_list)
print('Images: ', images_list)
print('Images count: ', images_count)

# Calculate the grid size:
grid_size = math.ceil(math.sqrt(images_count))

# Create plt plot:
fig, axes = plt.subplots(grid_size, grid_size, figsize=(result_figsize_resolution, result_figsize_resolution))

current_file_number = 0
for image_filename in images_list:
    x_position = current_file_number % grid_size
    y_position = current_file_number // grid_size

    plt_image = plt.imread(images_dir + '/' + images_list[current_file_number])
    axes[x_position, y_position].imshow(plt_image)
    print((current_file_number + 1), '/', images_count, ': ', image_filename)

    current_file_number += 1

plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
plt.savefig(result_grid_filename)