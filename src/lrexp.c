/* lrexp.c
 * computes the coupling term 
 *
 *   y_i = sum_{j=1:N} w_ij * f(x_j - x_i)
 *
 * in the most efficient way
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_sf.h>

#include "lrexp.h"
#include "datastruct.h"
#include "macros.h"

static const long TC[51][51] = { {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -3, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -8, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 5, 0, -20, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 18, 0, -48, 0, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -7, 0, 56, 0, -112, 0, 64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -32, 0, 160, 0, -256, 0, 128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 9, 0, -120, 0, 432, 0, -576, 0, 256, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 50, 0, -400, 0, 1120, 0, -1280, 0, 512, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -11, 0, 220, 0, -1232, 0, 2816, 0, -2816, 0, 1024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -72, 0, 840, 0, -3584, 0, 6912, 0, -6144, 0, 2048, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 13, 0, -364, 0, 2912, 0, -9984, 0, 16640, 0, -13312, 0, 4096, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 98, 0, -1568, 0, 9408, 0, -26880, 0, 39424, 0, -28672, 0, 8192, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -15, 0, 560, 0, -6048, 0, 28800, 0, -70400, 0, 92160, 0, -61440, 0, 16384, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -128, 0, 2688, 0, -21504, 0, 84480, 0, -180224, 0, 212992, 0, -131072, 0, 32768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 17, 0, -816, 0, 11424, 0, -71808, 0, 239360, 0, -452608, 0, 487424, 0, -278528, 0, 65536, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 162, 0, -4320, 0, 44352, 0, -228096, 0, 658944, 0, -1118208, 0, 1105920, 0, -589824, 0, 131072, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -19, 0, 1140, 0, -20064, 0, 160512, 0, -695552, 0, 1770496, 0, -2723840, 0, 2490368, 0, -1245184, 0, 262144, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -200, 0, 6600, 0, -84480, 0, 549120, 0, -2050048, 0, 4659200, 0, -6553600, 0, 5570560, 0, -2621440, 0, 524288, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 21, 0, -1540, 0, 33264, 0, -329472, 0, 1793792, 0, -5870592, 0, 12042240, 0, -15597568, 0, 12386304, 0, -5505024, 0, 1048576, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 242, 0, -9680, 0, 151008, 0, -1208064, 0, 5637632, 0, -16400384, 0, 30638080, 0, -36765696, 0, 27394048, 0, -11534336, 0, 2097152, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -23, 0, 2024, 0, -52624, 0, 631488, 0, -4209920, 0, 17145856, 0, -44843008, 0, 76873728, 0, -85917696, 0, 60293120, 0, -24117248, 0, 4194304, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -288, 0, 13728, 0, -256256, 0, 2471040, 0, -14057472, 0, 50692096, 0, -120324096, 0, 190513152, 0, -199229440, 0, 132120576, 0, -50331648, 0, 8388608, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 25, 0, -2600, 0, 80080, 0, -1144000, 0, 9152000, 0, -45260800, 0, 146227200, 0, -317521920, 0, 466944000, 0, -458752000, 0, 288358400, 0, -104857600, 0, 16777216, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 338, 0, -18928, 0, 416416, 0, -4759040, 0, 32361472, 0, -141213696, 0, 412778496, 0, -825556992, 0, 1133117440, 0, -1049624576, 0, 627048448, 0, -218103808, 0, 33554432, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -27, 0, 3276, 0, -117936, 0, 1976832, 0, -18670080, 0, 109983744, 0, -428654592, 0, 1143078912, 0, -2118057984, 0, 2724986880, 0, -2387607552, 0, 1358954496, 0, -452984832, 0, 67108864, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -392, 0, 25480, 0, -652288, 0, 8712704, 0, -69701632, 0, 361181184, 0, -1270087680, 0, 3111714816, 0, -5369233408, 0, 6499598336, 0, -5402263552, 0, 2936012800, 0, -939524096, 0, 134217728, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 29, 0, -4060, 0, 168896, 0, -3281408, 0, 36095488, 0, -249387008, 0, 1151016960, 0, -3683254272, 0, 8341487616, 0, -13463453696, 0, 15386804224, 0, -12163481600, 0, 6325010432, 0, -1946157056, 0, 268435456, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 450, 0, -33600, 0, 990080, 0, -15275520, 0, 141892608, 0, -859955200, 0, 3572121600, 0, -10478223360, 0, 22052208640, 0, -33426505728, 0, 36175872000, 0, -27262976000, 0, 13589544960, 0, -4026531840, 0, 536870912, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -31, 0, 4960, 0, -236096, 0, 5261568, 0, -66646528, 0, 533172224, 0, -2870927360, 0, 10827497472, 0, -29297934336, 0, 57567870976, 0, -82239815680, 0, 84515225600, 0, -60850962432, 0, 29125246976, 0, -8321499136, 0, 1073741824, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -512, 0, 43520, 0, -1462272, 0, 25798656, 0, -275185664, 0, 1926299648, 0, -9313976320, 0, 32133218304, 0, -80648077312, 0, 148562247680, 0, -200655503360, 0, 196293427200, 0, -135291469824, 0, 62277025792, 0, -17179869184, 0, 2147483648, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 33, 0, -5984, 0, 323136, 0, -8186112, 0, 118243840, 0, -1083543552, 0, 6723526656, 0, -29455450112, 0, 93564370944, 0, -218864025600, 0, 379364311040, 0, -485826232320, 0, 453437816832, 0, -299708186624, 0, 132875550720, 0, -35433480192, 0, 4294967296, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 578, 0, -55488, 0, 2108544, 0, -42170880, 0, 511673344, 0, -4093386752, 0, 22761029632, 0, -91044118528, 0, 267776819200, 0, -586290298880, 0, 959384125440, 0, -1167945891840, 0, 1042167103488, 0, -661693399040, 0, 282930970624, 0, -73014444032, 0, 8589934592, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -35, 0, 7140, 0, -434112, 0, 12403200, 0, -202585600, 0, 2106890240, 0, -14910300160, 0, 74977509376, 0, -275652608000, 0, 754417664000, 0, -1551944908800, 0, 2404594483200, 0, -2789329600512, 0, 2384042393600, 0, -1456262348800, 0, 601295421440, 0, -150323855360, 0, 17179869184, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -648, 0, 69768, 0, -2976768, 0, 66977280, 0, -916844544, 0, 8307167232, 0, -52581629952, 0, 240999137280, 0, -819082035200, 0, 2095125626880, 0, -4063273943040, 0, 5977134858240, 0, -6620826304512, 0, 5429778186240, 0, -3195455668224, 0, 1275605286912, 0, -309237645312, 0, 34359738368, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 37, 0, -8436, 0, 573648, 0, -18356736, 0, 336540160, 0, -3940579328, 0, 31524634624, 0, -180140769280, 0, 757650882560, 0, -2392581734400, 0, 5742196162560, 0, -10531142369280, 0, 14743599316992, 0, -15625695002624, 0, 12315818721280, 0, -6992206757888, 0, 2701534429184, 0, -635655159808, 0, 68719476736, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 722, 0, -86640, 0, 4124064, 0, -103690752, 0, 1589924864, 0, -16188325888, 0, 115630899200, 0, -601280675840, 0, 2334383800320, 0, -6880289095680, 0, 15547666268160, 0, -27039419596800, 0, 36108024938496, 0, -36681168191488, 0, 27827093110784, 0, -15260018802688, 0, 5712306503680, 0, -1305670057984, 0, 137438953472, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -39, 0, 9880, 0, -746928, 0, 26604864, 0, -543921664, 0, 7120429056, 0, -63901286400, 0, 411402567680, 0, -1960212234240, 0, 7061349335040, 0, -19502774353920, 0, 41626474905600, 0, -68822438510592, 0, 87841744879616, 0, -85678155104256, 0, 62646392979456, 0, -33221572034560, 0, 12060268167168, 0, -2680059592704, 0, 274877906944, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -800, 0, 106400, 0, -5617920, 0, 156900480, 0, -2677768192, 0, 30429184000, 0, -243433472000, 0, 1424085811200, 0, -6254808268800, 0, 21002987765760, 0, -54553214976000, 0, 110292369408000, 0, -173752901959680, 0, 212364657950720, 0, -199183403319296, 0, 140552804761600, 0, -72155450572800, 0, 25426206392320, 0, -5497558138880, 0, 549755813888, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, 41, 0, -11480, 0, 959728, 0, -37840704, 0, 857722624, 0, -12475965440, 0, 124759654400, 0, -898269511680, 0, 4808383856640, 0, -19570965872640, 0, 61508749885440, 0, -150732904857600, 0, 289407177326592, 0, -435347548798976, 0, 510407471005696, 0, -461013199618048, 0, 314327181557760, 0, -156371169312768, 0, 53532472377344, 0, -11269994184704, 0, 1099511627776, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                          {-1, 0, 882, 0, -129360, 0, 7537376, 0, -232581888, 0, 4393213440, 0, -55381114880, 0, 492952780800, 0, -3220624834560, 0, 15871575982080, 0, -60144919511040, 0, 177570714746880, 0, -411758179123200, 0, 752567256612864, 0, -1083059755548672, 0, 1219998345330688, 0, -1062579203997696, 0, 700809813688320, 0, -338168545017856, 0, 112562502893568, 0, -23089744183296, 0, 2199023255552, 0, 0, 0, 0, 0, 0, 0, 0},
                          {0, -43, 0, 13244, 0, -1218448, 0, 52915456, 0, -1322886400, 0, 21262392320, 0, -235521884160, 0, 1884175073280, 0, -11249633525760, 0, 51314117836800, 0, -181798588907520, 0, 505874334351360, 0, -1112923535572992, 0, 1940482062024704, 0, -2676526982103040, 0, 2901009890279424, 0, -2439485589553152, 0, 1557990796689408, 0, -729869562413056, 0, 236394999971840, 0, -47278999994368, 0, 4398046511104, 0, 0, 0, 0, 0, 0, 0},
                          {1, 0, -968, 0, 155848, 0, -9974272, 0, 338412800, 0, -7038986240, 0, 97905899520, 0, -963996549120, 0, 6988974981120, 0, -38370843033600, 0, 162773155184640, 0, -541167892561920, 0, 1423506847825920, 0, -2978414327758848, 0, 4964023879598080, 0, -6573052309536768, 0, 6864598984556544, 0, -5579780992794624, 0, 3454150138396672, 0, -1572301627719680, 0, 495879744126976, 0, -96757023244288, 0, 8796093022208, 0, 0, 0, 0, 0, 0},
                          {0, 45, 0, -15180, 0, 1530144, 0, -72864000, 0, 1999712000, 0, -35340364800, 0, 431333683200, 0, -3812168171520, 0, 25227583488000, 0, -128055803904000, 0, 507344899276800, 0, -1588210119475200, 0, 3959937231224832, 0, -7897310717542400, 0, 12604574741299200, 0, -16047114509352960, 0, 16168683558666240, 0, -12717552782278656, 0, 7638169839206400, 0, -3380998255411200, 0, 1039038488248320, 0, -197912092999680, 0, 17592186044416, 0, 0, 0, 0, 0},
                          {-1, 0, 1058, 0, -186208, 0, 13034560, 0, -484140800, 0, 11038410240, 0, -168586629120, 0, 1826663915520, 0, -14613311324160, 0, 88826010009600, 0, -418884762992640, 0, 1555857691115520, 0, -4599927086776320, 0, 10898288790208512, 0, -20758645314682880, 0, 31782201792135168, 0, -38958828003262464, 0, 37917148110127104, 0, -28889255702953984, 0, 16848641306132480, 0, -7257876254949376, 0, 2174833999740928, 0, -404620279021568, 0, 35184372088832, 0, 0, 0, 0},
                          {0, -47, 0, 17296, 0, -1902560, 0, 98933120, 0, -2967993600, 0, 57417185280, 0, -768506941440, 0, 7465496002560, 0, -54454206136320, 0, 305707823923200, 0, -1345114425262080, 0, 4699925501706240, 0, -13159791404777472, 0, 29693888297959424, 0, -54121865370664960, 0, 79611518093623296, 0, -94086339565191168, 0, 88551849002532864, 0, -65416681245114368, 0, 37078280867676160, 0, -15554790998147072, 0, 4547580092481536, 0, -826832744087552, 0, 70368744177664, 0, 0, 0},
                          {1, 0, -1152, 0, 220800, 0, -16839680, 0, 682007040, 0, -16974397440, 0, 283420999680, 0, -3363677798400, 0, 29544303329280, 0, -197734422282240, 0, 1030300410839040, 0, -4246086541639680, 0, 13999778090188800, 0, -37217871599763456, 0, 80146421910601728, 0, -140025932533465088, 0, 198181864190509056, 0, -226089827240509440, 0, 205992953708019712, 0, -147682003796361216, 0, 81414437990301696, 0, -33284415996035072, 0, 9499780463984640, 0, -1688849860263936, 0, 140737488355328, 0, 0},
                          {0, 49, 0, -19600, 0, 2344160, 0, -132612480, 0, 4332007680, 0, -91365980160, 0, 1335348940800, 0, -14192851599360, 0, 113542812794880, 0, -701176668487680, 0, 3405715246940160, 0, -13192098584985600, 0, 41159347585155072, 0, -104129631497486336, 0, 214414709191868416, 0, -359663383160553472, 0, 490450067946209280, 0, -540731503483551744, 0, 477402588661153792, 0, -332442288460398592, 0, 178383666978750464, 0, -71116412084551680, 0, 19826393672056832, 0, -3448068464705536, 0, 281474976710656, 0},
                          {-1, 0, 1250, 0, -260000, 0, 21528000, 0, -947232000, 0, 25638412800, 0, -466152960000, 0, 6034375680000, 0, -57930006528000, 0, 424820047872000, 0, -2432653747814400, 0, 11057517035520000, 0, -40383975260160000, 0, 119536566770073600, 0, -288405684905574400, 0, 568855350917201920, 0, -917508630511616000, 0, 1206989963132928000, 0, -1287455960675123200, 0, 1102487181118668800, 0, -746299014911098880, 0, 390051749953536000, 0, -151732604633088000, 0, 41341637204377600, 0, -7036874417766400, 0, 562949953421312},
                        };


/* takes state vector x and size N as input
 * and store meanx = mean(x) and range = range(x)
 */
int lrexpars(double *x, int N, double *meanx, double *range)
{
  int i;
	double minx = x[0],  maxx = x[0];
  *meanx = 0.0;
	for (i = 0; i < N; ++i) /* min/max x */
  {
    minx = minx > x[i] ? x[i] : minx;
    maxx = maxx < x[i] ? x[i] : maxx;
    *meanx += x[i];
  }
  *meanx /= (double)N;
  *range = (maxx - minx);

  /* fprintf(stderr,"meanx = %g, range = %g\n", *meanx, *range); */

  return 0;
}

/* lrexprank 
 * compute expansion rank 
 */
int lrexprank(coupling_function f,int N, int *p, double range)
{
  int         i;
  int         pmax = get_int("lrmax"); /* maximal polynomial order */ 
  double      sc = 0.0;
  int         oddity = 2; /* 0: even, 1: odd, 2: neither */
  int         pstep = 4; /* step on adaptive order p */
  double      abserr, 
              err = 0.0,
              pval,
              spe = 0.0;
  double      abstol = get_dou("lrabstol"),
              reltol = get_dou("lrreltol");
  size_t      errn = max(N/10,10);

  gsl_function     F;
  gsl_cheb_series *cs = gsl_cheb_alloc (pmax);

  F.function = f;
  F.params = &range;

  /* approximate the function f on [-1,1] */
  gsl_cheb_init (cs, &F, -1.0, 1.0);

  /* find out if f is even of odd */
  for (i = 0; i < pmax+1; i+=2)
    sc += fabs(cs->c[i]);

  if ( sc < (double)pmax/2 * 1e-8 ) /* this is almost odd function */
    oddity = 1;

  for (i = 1; i < pmax+1; i+=2)
    sc += fabs(cs->c[i]);

  if ( sc  < (double)pmax/2 * 1e-8 ) /* this is almost even function */
    oddity = 0;
  
  /* check whether function if even or odd */ 
  switch(oddity)
  {
    case 0 : /* start and stick with even rank */
      *p = 0;
      break; /* start and stick with odd rank */
    case 1 :  
      *p = -1;
      break;
    default : /* alternate even/odd rank */
      *p = -1;
      pstep = 3;
  }

  /* adaptative error on Chebychev approximation */
  do {
    err = 0.0;
    spe = 0.0;
    *p += pstep;
    for (i = 0; i < errn; ++i)
    {
      gsl_cheb_eval_n_err (cs, *p, -1.0+(double)i/(errn-1)*2.0, &pval, &abserr);
      err += abserr;
      spe += fabs(pval);
    }
  } while ( ( err/errn > (abstol + reltol*spe) ) & ( *p < pmax ) );

  gsl_cheb_free (cs);

  return 0;
}

/* compute_cheby_expansion 
 * compute Chebychev polynomial expansion of the
 * coupling function y[i] = 1/N*sum_j f(x[j] - x[i])
 * The algorithm is based on Chebychev approximation of the
 * function f(u) on u in [-range, range]
 * The algorithm runs in O(N*P^2), and is faster than the
 * direct method if N > P^2. Numerical tests show P up to 
 * 30. Advantageous if N > 2000
 */
int lrexp(coupling_function f, double *x, double *y, 
    int N, int p, double meanx, double range)
{
  int i,k,m,j,l;
  double phi;
  double *A= (double *)malloc( (p+1) * sizeof(double));
  double *B= (double *)malloc( (p+1) * sizeof(double));
  gsl_cheb_series *cs = gsl_cheb_alloc (p);

  gsl_function F;
  F.function = f;
  F.params = &range;

  gsl_cheb_init (cs, &F, -1.0, 1.0); /* approximate the function f on [-1,1] */
	double *cc = gsl_cheb_coeffs (cs); /* get the coefficients cc */
  /* The Chebychev series is computed using the convention
   *
   * f(x) = (c_0 / 2) + \sum_{n=1} c_n T_n(x)
   *
   * Since we want to sum
   *
   * c_0 + \sum_{n=1} c_n T_n(x)
   *
   * we need to set c_0 to c_0/2
   *
   */
  cc[0] /= 2;

  /* expand the Chebichev approximation  */
  for (m = p+1; m--;)  /* pre-compute the (p+1) moments A[m]: O(N*(p+1)) */ 
                       /* m goes from p to 0 */
	{	
		A[m] = 0.0;
		for (j = 0; j < N; ++j)  
		{
			A[m] += gsl_pow_int( (x[j] - meanx)/range, m);
		}
	}

  /* compute the weight c_k a_k,l */
  for (l = 0; l <= p; ++l)
  {
    B[l] = 0;
    for (k = l; k <= p; ++k)
    {
      B[l] += cc[k]*TC[k][l];
    }
  }
 
  for (i = 0; i < N; ++i) /* compute the coupling term: O(N*P^2)  */
	{
		y[i] = 0.0;
		for (m = 0; m <= p; ++m)
		{
			phi = 0.0;
      for (l = m; l <= p; ++l)
      {
        phi += B[l] * gsl_sf_choose ( l, m) * gsl_pow_int( -(x[i] - meanx)/range, l-m );
      }
			y[i] += A[m]*phi;
		}
    y[i] /= (double)N;
	}

  free(A);
  free(B);
	
  return 0;
}

/* lrexpw 
 * compute a polynomial approximation of the coupling term
 *
 *    y[i] = 1/N*sum_j W_ij f(x[j] - x[i])
 *
 * for i = 1,...,N
 *
 * where the NxN matrix W is of low rank r and is factorized 
 * 
 *    W = U*V
 *
 *    U Nxr matrix
 *    V rxN matrix
 * 
 *    W_ij = sum_m=1^r U_im V_mj
 * 
 * The algorithm is based on Chebychev approximation of the
 * function f(u) on u in [-range, range]
 * The algorithm runs in O(N*P^2), and is faster than the
 * direct method if N > P^2. Numerical tests show P up to 
 * 30. Advantageous if N > 2000
 *
 * the coupling term is return in the array y.
 *
 */
int lrexpw(coupling_function f, const double *U, const double *V, int r, double *x, double *y, 
    int N, int p, double meanx, double range)
{
  int i,k,m,j,l;
  int p1 = p+1;
  double phi;
  double *A       = (double *)malloc( (p1) * r * sizeof(double)); /* p+1 x r */
  double *B       = (double *)malloc( (p1) * sizeof(double));
  double ua;
  gsl_cheb_series *cs = gsl_cheb_alloc (p);

  gsl_function F;
  F.function = f;
  F.params = &range;

  gsl_cheb_init (cs, &F, -1.0, 1.0); /* approximate the function f on [-1,1] */
	double *cc = gsl_cheb_coeffs (cs); /* get the coefficients cc */
  /* The Chebychev series is computed using the convention
   *
   * f(x) = (c_0 / 2) + \sum_{n=1} c_n T_n(x)
   *
   * Since we want to sum
   *
   * c_0 + \sum_{n=1} c_n T_n(x)
   *
   * we need to set c_0 to c_0/2
   *
   */
  cc[0] /= 2;

  /* Moments A_km = sum_j=1^N V_mj psi_k(j)
   * k = 0...p
   * m = 1...r
   */
  for (k = 0; k < p1; ++k) 
	{	
    for (m = 0; m < r; ++m)
    {
		  *(A + p1*m + k) = 0.0;
		  for (j = 0; j < N; ++j)  
		  {
			  *(A + p1*m + k) += *(V + m*N + j)*gsl_pow_int( (x[j] - meanx)/range, k);
		  }
    }
	}

  /* B_l = sum_k=l^p c_k TC_kl 
   * l = 0,...,p
   */
  for (l = 0; l < p1; ++l)
  {
    B[l] = 0;
    for (k = l; k < p1; ++k)
    {
      B[l] += cc[k]*TC[k][l];
    }
  }
 
  for (i = 0; i < N; ++i) /* compute the coupling term: O(N*P^2)  */
	{
		y[i] = 0.0;
		for (k = 0; k < p1; ++k)
		{
      ua = 0;
      for (m = 0; m < r; ++m)
      {
        ua += *(U + i*r + m) * (*(A + p1*m + k));
      }
      phi = 0.0;
      for (l = k; l < p1; ++l)
      {
        phi += B[l] * gsl_sf_choose ( l, k) * gsl_pow_int( -(x[i] - meanx)/range, l-k );
      }
      y[i] += phi * ua;
		}
    y[i] /= (double)N;
	}

  free(A);
  free(B);
	
  return 0;
}

/* lrexpwp
 * low rank expansion for population-based system
 * compute a polynomial approximation of the coupling term
 *
 *    y[i] = 1/N*sum_j W_ij f(x[j] - x[i])
 *
 * for i = 1,...,N
 *
 * The algorithm is based on Chebychev approximation of the
 * function f(u) on u in [-range, range]
 * The algorithm runs in O(N*P^2), and is faster than the
 * direct method if N > P^2. Numerical tests show P up to 
 * 30. Advantageous if N > 2000
 *
 * the coupling term is return in the array y.
 *
 */
int lrexpwp(int iu, int iv, int r, coupling_function f, double *x, double *y, 
    int N, int p, double meanx, double range)
{
  int i,k,m,j,l;
  int p1 = p+1;
  double phi;
  double *A       = (double *)malloc( (p1) * r * sizeof(double)); /* p+1 x r */
  double *B       = (double *)malloc( (p1) * sizeof(double));
  double *U       = (double *)malloc( r * N * sizeof(double)); 
  double *Vs      = (double *)malloc( r * N * sizeof(double));
  double ua;
  par *myself_ = SIM->pop->start; 
  gsl_cheb_series *cs = gsl_cheb_alloc (p);

  /* construction of the matrices U and V */
  i = 0;
  while ( myself_ != NULL )
  {
    /* memcpy(void *restrict dst, const void *restrict src, size_t n); */
    memcpy(U + i, myself_->expr + iu, r * sizeof(double)); /* fill U row by row */
    memcpy(Vs + i, myself_->expr + iv, r * sizeof(double)); /* fill Vs row by row */
    myself_ = myself_->nextel;
    i += r;
  }

  gsl_function F;
  F.function = f;
  F.params = &range;

  gsl_cheb_init (cs, &F, -1.0, 1.0); /* approximate the function f on [-1,1] */
	double *cc = gsl_cheb_coeffs (cs); /* get the coefficients cc */
  /* The Chebychev series is computed using the convention
   *
   * f(x) = (c_0 / 2) + \sum_{n=1} c_n T_n(x)
   *
   * Since we want to sum
   *
   * c_0 + \sum_{n=1} c_n T_n(x)
   *
   * we need to set c_0 to c_0/2
   *
   */
  cc[0] /= 2;

  /* Moments A_km = sum_j=1^N V_mj psi_k(j)
   * k = 0...p
   * m = 1...r
   */
  for (k = 0; k < p1; ++k) 
	{	
    for (m = 0; m < r; ++m)
    {
		  *(A + p1*m + k) = 0.0;
		  for (j = 0; j < N; ++j)  
		  {
			  *(A + p1*m + k) += *(Vs + m + j*r)*gsl_pow_int( (x[j] - meanx)/range, k);
		  }
    }
	}

  /* B_l = sum_k=l^p c_k TC_kl 
   * l = 0,...,p
   */
  for (l = 0; l < p1; ++l)
  {
    B[l] = 0;
    for (k = l; k < p1; ++k)
    {
      B[l] += cc[k]*TC[k][l];
    }
  }
 
  for (i = 0; i < N; ++i) /* compute the coupling term: O(N*P^2)  */
	{
		y[i] = 0.0;
		for (k = 0; k < p1; ++k)
		{
      ua = 0;
      for (m = 0; m < r; ++m)
      {
        ua += *(U + i*r + m) * (*(A + p1*m + k));
      }
      phi = 0.0;
      for (l = k; l < p1; ++l)
      {
        phi += B[l] * gsl_sf_choose ( l, k) * gsl_pow_int( -(x[i] - meanx)/range, l-k );
      }
      y[i] += phi * ua;
		}
    y[i] /= (double)N;
	}

  free(A);
  free(B);
  free(U);
  free(Vs);
	
  return 0;
}
/* lrkern
 * **Main user function**
 * computes the coupling term 
 *
 *   y_i = sum_{j=1:N} f(x_j - x_i)
 *
 * in the most efficient way
 */
int lrkern(coupling_function f, double *x, double *y, int N)
{
  int         p; /* polynomial order */ 
	double      range,
              meanx;
  
  lrexpars(x, N, &meanx, &range); /* get meanx, range */
  lrexprank(f,N,&p, range); /* get rank p */
  lrexp(f,x,y,N,p,meanx,range); /* compute high accuracy */

  return 0;
}

/* lrwkern
 * **Main user function**
 * computes the coupling term 
 *
 *   y_i = sum_{j=1:N} W_ij f(x_j - x_i)
 *
 * where W_ij can be factorized as U_i*V_j
 * in the most efficient way
 */
int lrwkern(const double *U, const double *V, int r, coupling_function f, double *x, double *y, int N)
{
  int         p; /* polynomial order */ 
	double      range,
              meanx;
  
  lrexpars(x, N, &meanx, &range); /* get meanx, range */
  lrexprank(f,N,&p, range); /* get rank p */
  lrexpw(f,U,V,r,x,y,N,p,meanx,range); /* compute high accuracy */

  return 0;
}

/* lrwpkern
 * **Main user function**
 * computes the coupling term 
 *
 *   y_i = sum_{j=1:N} W_ij f(x_j - x_i)
 *
 * where W_ij can be factorized as U_i*V_j
 * in the most efficient way
 */
int lrwpkern(int iu, int iv, int r, coupling_function f, double *x, double *y, int N)
{
  int         p; /* polynomial order */ 
	double      range,
              meanx;
  
  lrexpars(x, N, &meanx, &range); /* get meanx, range */
  lrexprank(f,N,&p, range); /* get rank p */
  lrexpwp(iu,iv,r,f,x,y,N,p,meanx,range); /* compute high accuracy */

  return 0;
}

