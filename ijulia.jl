using Base64
using PyPlot
using PyCall
animation = pyimport("matplotlib.animation")

# Setup fonts and styles
rcParams = PyDict(matplotlib["rcParams"])
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = "Routed Gothic"
rcParams["font.size"] = 15
rcParams["axes.linewidth"] = 1
rcParams["scatter.marker"] = "o"
rcParams["xtick.direction"] = "in"
rcParams["ytick.direction"] = "in"

function animate(f, frames, figsize, interval=1000/25.)
  filename = tempname() * ".mp4"
  fig = PyPlot.figure(figsize=figsize)
  withfig(fig) do
    myanim = animation.FuncAnimation(fig, f, frames=frames, interval=interval)
    myanim[:save](filename, bitrate=-1, extra_args=["-vcodec", "libx264", "-pix_fmt", "yuv420p"])
  end
  data = base64encode(open(filename))
  display("text/html", """<video controls src="data:video/x-m4v;base64,$data">""")
end

nothing