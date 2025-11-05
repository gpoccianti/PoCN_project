# PoCN_project

In this repository, you will find reports, code, and data output for the projects related to the course *Physics of Complex Networks: Structure and Dynamics* (AY 2023/2024, Master degree program in *Physics of Data*, University of Padova).

## Chosen Projects:

### 1. Theoretical
- **20:** Swarmalators (different phases) [score: 0.6]
  > "Swarmalators" are theoretical entities, first introduced by O'Keeffe, Hong, and Strogatz (2017), that bridge the gap between synchronization and swarming phenomena. They are defined as oscillators whose internal phase dynamics and spatial movements are bidirectionally coupled: their phase affects how they move, and their position affects how they synchronize.
  > 
  > The objective of this project is to implement and simulate the foundational swarmalator model. We numerically analyze its behavior by exploring the key **K-J parameter space** to reproduce the five distinct collective states identified in the literature: **Static Sync**, **Static Async**, **Static Phase Wave**, **Splintered Phase Wave**, and **Active Phase Wave**. The analysis also investigates the model's dynamics across different synthetic network topologies, as required by the course.

<p align="center">
  <a href="/data/task_2Z0/Animations/Complete Network/swarmalators_complete_N=300,T=2000,dt=0.05,K=-0.60,J=0.90.mp4" title="Click to watch the full video">
    <img src="/data/task_20/Animations/Complete Network/swarmalators_complete_N=300,T=2000,dt=0.05,K=-0.60,J=0.90.gif" alt="Swarmalators Animation Preview">
  </a>
  <br>
  <em>N=300,T=2000,dt=0.05,K=-0.60,J=0.90</em>
</p>

- **29:** Language competition dynamics [score: 0.4]
  > This project simulates a well-known dynamical model from sociophysics, based on the paper "Modelling the dynamics of language death" by Abrams & Strogatz (2003).
  >
  > The model describes the competition between two languages (X and Y) competing for speakers. The core assumption is that the attractiveness of a language increases with both its number of speakers and its perceived "status" ($s$), which reflects social or economic opportunities.
  >
  > The dynamics are governed by the equation $\frac{dx}{dt} = yP_{yx}(x,s) - xP_{xy}(x,s)$, where $x$ is the fraction of speakers of language X. In its basic form, the model predicts that two languages cannot coexist stably; one will inevitably drive the other to extinction.
  >
  > The objective is to computationally simulate this model, analyze the influence of its key dynamical parameters (especially the **status $s$**), and study how these dynamics behave across various synthetic network topologies, as required by the course.

### 2. Data Analytics
- **43:** Social connectedness index I from Facebook [score: 1.0]
  > This project is a data analytics task focusing on the **Facebook Social Connectedness Index (SCI)**, a dataset that quantifies the relative probability of friendship ties between two geographic areas.
  >
  > The primary objective is to construct, visualize, and analyze the internal social network of individual **US states**. For each state, a network is built where:
  > * **Nodes** are the counties within that state, with attributes for `nodeID`, `nodeLabel` (county name), `latitude`, and `longitude`.
  > * **Edges** represent the SCI "friendship" link between pairs of counties (`nodeID_from`, `nodeID_to`).
  >
  > After constructing the graph for each state, the project involves a **simplified comparative network analysis**. This includes calculating and comparing key topological metrics (e.g., density, average degree, clustering) to understand and contrast the structural patterns of social connectedness across different US states.