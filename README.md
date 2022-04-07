# L00_DosDim123_lenFn
This code is uploaded by authorization and request of the first author and corresponding author.


## Introduction

Pulmonary fissures are double layers of invaginations of visceral pleura that anatomically separate the lungs into lobes and segments. Interlobar fissures are the physical boundaries between pulmonary lobes, dividing the human lung into five lobes (three in the right lung and two in the left lung). Accessory fissures often occur between bronchopulmonary segments but may also enter subsegmental or interbronchial planes. Pulmonary fissures are important landmarks for delineation of pulmonary anatomy, and have significant value in localization of lesions and assessment of disease processes. Particularly, there is an increasing requirement for quantification of fissure integrity or completeness, which is closely relevant to lung disease characterization and can be calculated as the percentage of the interlobar border defined by a fissure. However, automated or computer-aided segmentation of pulmonary fissures in CT images is not an easy task. The main challenges come from complicating factors such as their thin, weak and variable structures, pathological deformation, inhomogeneous intensity, imaging noise along with interferences from adjacent vessels, bronchi and pathological structures (e.g. fibrotic tissue).

Several automatic or semi-automatic methods have been reported in the literature to realize pulmonary fissures segmentation. With the lung parenchyma segmentation, pulmonary airway and vessel trees extraction techniques gradually becoming mature, an indirect strategy is to utilize these detected tissues to guess the fissure position on the basis of some prior assumptions. Among them, the sparse distribution of vessels and bronchi in the neighborhood of fissures is commonly used anatomical knowledge, which has been integrated under different segmentation frameworks like the watershed transform, Voronoi division , adaptive sweeping , minimal path  and neural network classifier . Another more global constraint is integrated in lung atlas schemes including the single atlas search initialization and a multi-atlas selection mechanism recently proposed by van Rikxoort et al. Although these indirect methods could provide an approximate estimation of the fissure position and thus help to reduce computational burden, the most reliable information still comes from the object itself. This is also the reason why the fissure appearance and shape characteristics are attracting more attention for accurate localization and refinement. 

Much effort has been made to efficiently exploit the direct fissure information. Based on the fact that the profile of pulmonary fissures across transverse planes can be approximated with piecewise straight lines, Kubo et al.  presented a 2D VanderBrug operator to enhance the linear structures and simultaneously suppress streak artifacts. Their method was later extended to 3D space as a sheet-emphasis filter. Zhang et al.  began with a 2D ridgeness measure to enhance fissure contrast, the direction and intensity continuity was integrated under a fuzzy reasoning framework to sift out the fissure objects. Wiemker et al.  proposed two equivalent fissure likelihood filters using structure tensors and Hessian matrices for shape description. A similar likelihood function was recently defined by Lassen et al.  and the plane filter of Li et al.  can also be considered a pilot method to enhance the fissure-like tissues. To automatically distinguish what should be enhanced or suppressed, van Rikxoort et al.  presented a supervised filter, which adopts a pattern recognition technique to select the most important features and was further harnessed to group the unordered structures into bigger fissure plates . A computational geometry approach has been developed by Pu et al. where they utilized a statistical approach to extract the pulmonary fissures in 3D space after iteratively smoothing the triangle meshes with a Laplacian filter. The same authors later suggested an anisotropic morphological operator as a post-processing method to smoothen fissure surfaces and fill small holes. Furthermore, their method was improved by a piecewise plane fitting algorithm to directly identify fissure patches from the original lung sub-volumes. Ross et al.  employed a particle system that samples the image domain and provides a set of candidate fissure locations based on image derivative features. A maximum a posteriori estimation was then applied to eliminate poor candidates and remove residual noise particles. Later, their ridge surface sampling scheme was merged with a lobe boundary shape model to improve fissure discrimination. Recently, Kinder et al.  proposed a line enhancing filter as an extension to the previous Hessian filter. Their method is close to our proposed filter but they used a single stick template. 

Motivated by a line detection model in speckle images, we present a derivative of stick (DoS) filter for fissure enhancement with emphasis on direction estimation and interference suppression. The basic idea is to probe the presence of 2D fissure profiles across section planes by formulating it as a M-ary hypothesis testing problem, where a rectangle neighborhood is divided into various straight line-segments (i.e. sticks) corresponding to potential fissure directions. Different from traditional isotropic or spot kernel filters, our method defines a likelihood measure by investigating appearance feature along orienting templates, which make it more appropriate for extremely anisotropic and thin elongated structures. For subsequent segmentation, we introduce a postprocessing pipeline based on connected component analysis, where multi-threshold binarized images are directly merged after removing adhering clutters with a branch-point detection algorithm. 

In this paper, our main purpose was to extract pure fissure patches rather than generate a full lobe segmentation. This means that only the visible fissures on CT are extracted and no interpolation operation was applied to extend the fissure plane or fill its inner holes. An early version of this method has been presented at a conference . In this current work, the algorithms are improved further and the experiments have been extended into a full validation. The remainder of the paper is organized as follows. In Section II, the methods are described in details. The data and reference are described in Section III. We give the experiment and evaluation results in Section IV, and in Section V the conclusions are presented. 
