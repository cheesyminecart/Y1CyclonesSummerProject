import imageio.v2 as imageio
import os

png_folder = r"C:\Users\Jack\OneDrive - Imperial College London\animation frames isabel" #This is the folder with all the frames in it
output = r"C:\Users\Jack\OneDrive - Imperial College London\isabel.mp4" #Where you want to save the animation

frames = []

for i in range(1, 55): #change 55 to the number of frames + 1 (e.g. if there are 34 frames the RHS number should be 35)
    filename = os.path.join(png_folder, f"frame_{i:02d}.png")
    frames.append(imageio.imread(filename))

imageio.mimsave(output, frames, fps=10)
