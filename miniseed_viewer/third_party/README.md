# Bundled third-party code

All libraries are copied in as source so the build never needs the network
and never changes underneath us.

| Library | Version | Use | License |
|---|---|---|---|
| libmseed | 3.0.0 | reads miniSEED | Apache 2.0 |
| Dear ImGui | 1.91.9b | user interface | MIT |
| ImPlot | commit 851cdbd | plots | MIT |
| GLFW | 3.4 | window + mouse/keyboard | zlib |
| pocketfft (header-only C++) | commit c90e55b | FFT | BSD-3 |
| Roboto-Medium.ttf (embedded in fonts/roboto_medium.h) | from ImGui 1.91.9b misc/fonts | UI font | Apache 2.0 |
