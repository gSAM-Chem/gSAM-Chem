ffmpeg -framerate 10 -i ./image.%06d.png -c:v libx264 -pix_fmt yuv420p output.mp4
