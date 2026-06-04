ffmpeg -framerate 20 -i ./image.%06d.png -c:v libx264 -pix_fmt yuv420p output.mp4
