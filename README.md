# Joint-User-Association-and-In-band-Backhaul-Scheduling-and-in-5G-mmWave-Networks
Matlab Simulations for 

"Joint Load Balancing and Interference Mitigation in 5G Heterogeneous Networks," IEEE Transactions on Wireless Communications 16 (9), 6032 - 6046, Sept. 2017.

URL: http://ieeexplore.ieee.org/abstract/document/7959611/

"Joint In-Band Backhauling and Interference Mitigation in 5G Heterogeneous Networks," European Wireless 2016; 22th European Wireless Conference, Oulu, Finland, 2016, pp. 1-6. 

URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&amp;arnumber=7499273&amp;isnumber=7499250



The package contains a outline simulation, based on Matlab, that helps reproduce all the numerical results and figures in the article. We encourage you to also perform reproducible research!

# Abstract I
We study the problem of joint load balancing (user association and user scheduling) and interference mitigation (beamforming design and power allocation) in 5G heterogeneous networks (HetNets) in which massive multiple-input multiple-output (MIMO) macro cell base station (BS) equipped with a large number of antennas, overlaid with wireless self-backhauled small cells (SCs) are assumed. Self-backhauled SC BSs with full-duplex communication employing regular antenna arrays serve both macro users and SC users by using the wireless backhaul from macro BS in the same frequency band. We formulate the joint load balancing and interference mitigation problem as a network utility maximization subject to wireless backhaul constraints. Subsequently, leveraging the framework of stochastic optimization, the problem is decoupled into dynamic scheduling of macro cell users, backhaul provisioning of SCs, and offloading macro cell users to SCs as a function of interference and backhaul links. Via numerical results, we show the performance gains of our proposed framework under the impact of small cells density, number of base station antennas, and transmit power levels at low and high frequency bands. We further provide insights into the performance analysis and convergence of the proposed framework. The numerical results show that the proposed user association algorithm outperforms other baselines. Interestingly, we find that even at lower frequency band the performance of open access small cell is close to that of closed access at some operating points, the open access full- duplex small cell still yields higher gain as compared to the closed access at higher frequency bands. With increasing the small cell density or the wireless backhaul quality, the open access full- duplex small cells outperform and achieve a 5.6x gain in terms of cell-edge performance as compared to the closed access ones in ultra-dense networks with 350 small cell base stations per km2.

# Abstract II
In this paper, we study the problem of joint in-band backhauling and interference mitigation in 5G heterogeneous networks (HetNets) in which a massive multiple-input multiple-output (MIMO) macro cell base station equipped with a large number of antennas, overlaid with self-backhauled small cells is assumed. This problem is cast as a network utility maximization subject to wireless backhaul constraints. Due to the non-tractability of the problem, we first resort to random matrix theory to get a closed-form expression of the achievable rate and transmit power in the asymptotic regime, i.e., as the number of antennas and users grows large. Subsequently, leveraging the framework of stochastic optimization, the problem is decoupled into dynamic scheduling of macro cell users and backhaul provisioning of small cells as a function of interference and backhaul links. Via simulations, we evaluate the performance gains of our proposed framework under different network architectures and low/high frequency bands. Our proposed HetNet method achieves the achievable average UE throughput of 1.7 Gbps as well as ensures 1 Gbps cell-edge UE throughput when serving 200 UEs per km2 at 28 GHz with 1 GHz bandwidth. In ultra-dense network, the UE throughput at 28 GHz achieves 62x gain as compared to 2.4 GHz.

# System requirement*
-  MATLAB version: 9.1 (R2016b)  https://www.mathworks.com/downloads/

- OS: Windows 7 amd64 version 6.1

- Java version: 1.7.0_60

- Yalmip lastest version https://yalmip.github.io/?n=Main.Download

- Mosek solver version 7 at least : https://www.mosek.com/

# Acknowledgements
The authors would like to thank the Finnish Funding Agency for Technology and Innovation (Tekes), Nokia, Huawei, and Anite for project funding. The Academy of Finland funding through the grant 284704 is also acknowledged.

# License and Referencing
This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
