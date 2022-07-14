# xkmeans
This is a pure C++ implementation about the traditional k-means and a dozens of its variants.

The performance evaluation and the details of our proposals of k-means variants, such as k-means#, k-sums, and sequential-k-sums can be found from our papers.
###### [1] Wan-Lei Zhao, Cheng-Hao Deng, Chong-Wah Ngo: k-means: a revisit, Neurocomputing, 2018.
###### [2] Wan-Lei Zhao, SY Lan, Run-Qing Chen, Chong-Wah Ngo: k-sums clustering: a stochastic optimization approach, CIKM 2021.
###### [3] Cheng-Hao Deng, Wan-Lei Zhao: Fast k-Means Based on k-NN Graph, ICDE 2018.


# The k-means variants
##### k-means
##### k-means++
##### mini-Batch
##### hartigan method
##### k-means^# (our propoal, see paper [1])
##### graph-based k-means (our proposal, see paper [3])
##### Repeated bisecting k-means (from Ying Zhao and Geoge Karypis [4])
##### sequential k-means
##### sequential k-sums (our proposal, see paper [2])
##### The bisecting version of above k-means variants have been also integrated.

###### [4] Ying Zhao, Geoge Karypis: Criterion functions for document clustering: Experiments and analysis, 2001.

# Compile and Install
#### This project is implemented with C++. It can be compiled smoothly under Ubuntu (MacOS) with g++ 8.0 or later. Go to under the directory, and run
### <b> make release </b>

#### In Windows, it can be smoothly compiled with MingW.

# Authors
###### Ranked by the contributed number of code lines
##### Wan-Lei Zhao (2010 - now)
##### Cheng-Hao Deng (2015 - 2017)
##### Run-Qing Chen (2018 - 2019)
##### Hui Ye (2018 - 2019)

Great thanks to my students Mr. Cheng-Hao Deng, Mr. Run-Qing Chen, and Miss Hui Ye!!! This project would not be in its current shape without any of you.

# History
The earliest codes of this project could be traced back to 2010, right after Wan-Lei Zhao received his Ph.D in City University of Hong Kong. He was a little bit free at that time. He tried to re-implement RB k-means (proposed Ying Zhao and George Karypis from UMD), which is very efficient to build visual vocabulary. At that time, SIFT feature was widely used for image search tasks. It is necessary to call k-means to build visual word vocabulary based on SIFT feature. RB k-means with no doubt is the best option for its high efficiency. So I re-implemented RB k-means at that time and the traditional k-means as well.

Wan-Lei Zhao joined with Xiamen University in Sep. 2014. Wan-Lei Zhao had my first gradudate student Cheng-Hao Deng in 2015. At that time, he wanted to find out an efficient way to build the hyerlinks between images. Cheng-Hao Deng made a breif survey, and discussed with Wan-Lei Zhao "whether it would be efficient if we use clustering approach". Then Wan-Lei Zhao went back to sort out the codes, and working on it together with Cheng-Hao Deng. Then this project was developped and grows day by day.
