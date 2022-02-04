# SSBI-effect

Over the past few years, we have witnessed the transfer of a huge amount of information through data centers, mainly due to the use of social networks and cloud computing. Already in 2016 90% of global internet traffic went through a central information center. In 2019, the central traffic between data reached a global capacity of Exabyte 905. Therefore, it is necessary to develop optical technologies for the transfer of information at a very high transmission rate of 100 Gb / s (or more) cheaply.

The aim of the project is to investigate the operations and performances of a direct detection optical receiver for coherent QAM optical modulation methods as a candidate for cheap and efficient technology for the purpose of high-bit rate data center interconnect. In the receiver after the sampler, we will implement the algorithm that cancels the SSBI iteratively for proper detection of QAM methods.

The world demands that technology knows how to transmit a huge amount of information in as short a time as possible - one way to implement this requirement is through optical communication.
In the project, we are dealing with an optical communication system at a high transmission rate of 60 Gbaud. We will refer to the ratio of optical fiber as chromatic dispersion, and solving these solution.
Next, we would like to study the performance of a direct receiver - a high-speed signal is transmitted with a carrier wave at one end of the signal width (like an analog modulation method - SSB), after the direct detection, we can get a non-linear called SSBI, which is actually an unwanted signal .
 
We will examine the results obtained and offer some solutions for better quality detection of the signal.

The project is the basic simulation in the Matlab environment, at each stage a performance investigation is made and conclusions are drawn regarding the next steps for the continuation of the project.
 
 
 ---
 

### **For exmaple: Task No. 4 - Direct Rx - Algorithm for SSBI cancelation , without using clipping:**


![image](https://user-images.githubusercontent.com/92098070/152528839-136ddb75-ba44-4478-8583-437bf103190d.png)


![image](https://user-images.githubusercontent.com/92098070/152529172-41fd06a2-cf8a-4e63-9097-10e3aca8ad2d.png)



It can be seen that the increase in attractions leads to instability of the system.
Later we will use clipping and offer solutions to this phenomenon.
