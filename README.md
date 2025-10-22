# Vertex-alignment-and-changepoint-localization-in-network-time-series

This repository contains code supporting the text "Vertex-alignment-and-changepoint-localization-in-network-time-series" by Tianyi Chen[^1], Mohammad Sharifi Kiasari[^2], Sijing Yu, Youngser Park, Avanti Athreya, Vince Lyzinski, Carey Priebe, Zachary Lubberts

## Abstract
Changepoint localization in a time series of networks generally relies on accurate vertex correspondence between
network realizations at different times. However, such vertex alignments are often misspecified or even unknown.
We show how misalignment between network realizations at different times can weaken their underlying correlation,
impeding inference procedures that rely on accurate correlation estimates. We explore the conceptual relationship
between graph matching and optimal transport, an approach often considered for mitigating errors from misalignment.
We construct two illustrative models for network evolution, each with a similar changepoint. We compare techniques for
changepoint localization, ranging from the simple network statistic of average degree to the more involved and recently
developed procedure of Euclidean mirrors. In one model, vertex misalignment causes comparatively little error, and
in the other, it seriously impairs localization, a problem where graph matching proves ineffective. Despite this, the
Euclidean mirror procedure can still extract meaningful signal when the misalignment portion is small. Finally We
present simulations to illustrate these contrasting effects on approaches to localization.


## Repository Structure

- **figures** — Includes code to reproduce all figures presented in the paper exactly.  
- **simulation** — Provides the full implementation and scripts used to generate the experimental results reported in the paper, multiple functions adopted from [TianyiChen97 / Euclidean-mirrors-and-first-order-changepoints-in-network-time-series](https://github.com/TianyiChen97/Euclidean-mirrors-and-first-order-changepoints-in-network-time-series)


[^1]: Co-first author
[^2]: Co-first author

## Citation

Chen, T., Lubberts, Z., Athreya, A., Park, Y., & Priebe, C. E. (2025).  
*Euclidean Mirrors and First-Order Changepoints in Network Time Series* [Computer software]. GitHub.  
Available at: [https://github.com/TianyiChen97/Euclidean-mirrors-and-first-order-changepoints-in-network-time-series](https://github.com/TianyiChen97/Euclidean-mirrors-and-first-order-changepoints-in-network-time-series)

