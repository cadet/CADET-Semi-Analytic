// =============================================================================
//  CADET-semi-analytic - The semi-analytic extension of CADET
//  
//  Copyright © 2015-present: Samuel Leweke¹² and the CADET-Semi-Analytic
//  authors, see the AUTHORS file
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//    ² University of Cologne, Cologne, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>

#include <vector>
#include <limits>
#include "MPReal.hpp"
#include "BesselZeros.hpp"

namespace
{
	inline bool equal(const mpfr::mpreal& ref, const mpfr::mpreal v)
	{
		return (abs(ref - v) <= std::numeric_limits<mpfr::mpreal>::epsilon()) || (abs(ref - v) / ref <= std::numeric_limits<mpfr::mpreal>::epsilon());
	}
}

TEST_CASE("Bessel zeros", "[Bessel],[CI]")
{
	mpfr::mpreal::set_default_prec(mpfr::digits2bits(100));

	// Reference values obtained by Mathematica (100 digits precision)
	const std::vector<mpfr::mpreal> ref
	{
		"0",
		"3.831705970207512315614435886308160766564545274287801928762298989918839309519011470214112874757423127",
		"7.015586669815618753537049981476524743276311502911313896055377826985496015502018663072714930179466458",
		"10.17346813506272207718571177677584406981951250019168555611465006811578704378288387382891893264510929",
		"13.32369193631422303239368412694787675121664473135786578547757152649656706334730478254711901794883995",
		"16.47063005087763281255246047098955144943812682227312576994442000794444256764074389900293149702940015",
		"19.61585851046824202112506588413750985024740266188054464735144476447080110172729827148538916115751770",
		"22.76008438059277189805300515218225759290537073807322687200507713025427356258966419457032670230020112",
		"25.90367208761838262549585544597987428790542703136724764136710449434505671416316505888147226836787887",
		"29.04682853491685506664781988353196110041417179308387566604012557337657839553619620936782533808676893",
		"32.18967991097440362662298410446036921905286771101497692036035465449720252920020026614512664228153970",
		"35.33230755008386510263447902251901446982934885987796525950522936960467593631437701003475084092397415",
		"38.47476623477161511205219755771661255497612689214613674271060785101617366150762879331639878430043424",
		"41.61709421281445088586351680506028920015120551602838853428546967173558795423254302975468298286226610",
		"44.75931899765282173277935271321214421818367022549847169807529030334648501424451866130043148826969685",
		"47.90146088718544712127400872250750718946553530195782853593514817141513730553135992482551538524868981",
		"51.04353518357150946873303463322406072863248537382081368603368240738292513466632428730738834243673259",
		"54.18555364106132053209996621453388893786434688422355844664222120864833562548095773069702109143015621",
		"57.32752543790101074509050424375052534918564187373742515421581109399737200548403490614447492297841187",
		"60.46945784534749155939874980838314984863581942099681168326949017495224472312920586186020239262221351",
		"63.61135669848123263103976241787362606002869218523134601188005677055284043560862035724798197495538646",
		"66.75322673409849341530525975004239661319910178809906444990721542666543857130449167050712063122109012",
		"69.89507183749577396973053643549981188025788331273032763768073657765686393275790203922206937960992295",
		"73.03689522557383482650611756909204894697544929328739503832140209107973677753472721548765138297033489",
		"76.17869958464145757285261462353460715231201091394723146404802976342382762150600887213641982467131467",
		"79.32048717547629939118448487248769143631435346876750254407144138836072227571884195128657049918183144",
		"82.46225991437355645398661064878121486353823164041695797183781769364355834725533646438979457489068860",
		"85.60401943635023096594942549338001602988691808389713000384015905516751370675059304829461054596430131",
		"88.74576714492630690373591643485379686825991667375800516990733292403446115583720441845274994963392240",
		"91.88750425169498528055362221449036745244349342099263804862879067639258167473359117131717385043360127",
		"95.02923180804469526805099818717409626011528711318811803767845554158691239567892359552113579828502756",
		"98.17095073079078197353775916085095547507022596313898004527267042854069499829903932785268654383992004",
		"101.3126618230387301371410563886518962183802038506392348094490992159908704032865312739837734408982717",
		"104.4543657912827600713634281396138254919369993599393141126066596193530268567390833895614890665369275",
		"107.5960632595091721826703642776082427270518710170091435863961143169149688760921353349602515795300942",
		"110.7377547808992151086086528882720795026222379407361112321239197453891134696598946723350550927143264",
		"113.8794408475949981348841749284283841537444709037317037034600551625807799280127121331591471220249397",
		"117.0211218988924250275764946014608412039486324166026679867125622122171325563177920727415641060596697",
		"120.1627983281490037581194078291668284720100287858882494598318814582493671041742054954934200801145477",
		"123.3044704886357180167600320687683589634812305979296099059542754415543098383205134604452527841013184",
		"126.4461386985165956977944804958388063069254245114358459781624291381733222765786006432972305427214918",
		"129.5878032451039967537414178413553395608656031039302054675356623627614241367468600697311874560977458",
		"132.7294643885096158867745973517548219009715509635320794070520678625168676248212438588987148372017829",
		"135.8711223647890005918015682194598132859015500176100124847826161607243591617589615011600945078592978",
		"139.0127773886597041784335461359553558207479502934947712140641156755844982700565653597089956999343913",
		"142.1544296558590290327009080997622490319967030604518025824791989211441613757662974947912739833774351",
		"145.2960793451959072324221508550143882477650791242003823672798171853306251329946994600483697959210646",
		"148.4377266203422303959392770262723049317877931274895843206145087742951216895976832088190068545183243",
		"151.5793716314014279927835042222342353964576147819375424732107216496779084915520225125461970178659982",
		"154.7210145162859535247665556518388510442927579133410356162762765327570831586654577022425102767524954",
		"157.8626554019302978050946696086597330486145694405456494313108221671998665135323890078877005655837640",
		"161.0042944053619934638934154090867336444933314748304785129740059367598650989580640215341705167961628",
		"164.1459316346496354021325267799765728809173970888529652836528898912646773677646226692518351393778403",
		"167.2875671897440838035648478578932261985025893145292612810973711441093578406368101413257383730814612",
		"170.4292011632266323477454974197510711799585043261854457753509983106032141573431201484922031614392788",
		"173.5708336409759286303670408628485014933839179115839098348903662915734967773940802753233443287266626",
		"176.7124647027637574552977373410500649046297593262796211072305931842162724258995493818517258484312698",
		"179.8540944227883848450849029063305509672004580314388248187730558148459261689771824640565797169605563",
		"182.9957228701529660840843730950345005378658428153534326851535590300987650420276371440723104134343974",
		"186.1373501092955080202398480948943831236625334277064463786412526661805053162277090951048623281922432",
		"189.2789762003760140932258006787760093269353412063299474494797083922206191961800136235526972808042478",
		"192.4206011996257054217711315547800983214950097494351189379069785731913693907806412723527255898562149",
		"195.5622251596625824307829356041230559925657061350128780785989473037017736439121132542105025372435343",
		"198.7038481297770521261181065789856227571481156562985992184465351965704977625497668332563927410219747",
		"201.8454701561908823049998174349850182993306599367085241557455228926507597120955223013891864457938989",
		"204.9870912822923441443587342300625607375308682940723998736780842850180458861154588610331595197995474",
		"208.1287115488500590814872465223503123340122307521193104863844542740448942957202225440814722458328406",
		"211.2703309942077666144627914099820645818711135187645248275416545480492580244204126740574882776823570",
		"214.4119496544619698287002708497179120215894014688857388524202694542417876957243371311014126034253564",
		"217.5535675636241894017847905098618069260469325521222435060828700908286351747119330559673361727903231",
		"220.6951847537693597448129102834219122803056234942289655674697585259215213725993670278417390584169959",
		"223.8368012551717287402915535098319723149254701554956333297249337091062426499989817532462681730401899",
		"226.9784170964294717884801879548191462793766329904583154939840153711676882602077336299060725500424178",
		"230.1200323045790986476247695568069438495838331396776745459195739785528738570019336720083939456725588",
		"233.2616469052006153543463062410181813647518202190889401184536785381801019165117011648838233789914575",
		"236.4032609225143012086734229512145946319448029124450282156128539189649824924724313138902593764525035",
		"239.5448743794698705812588571641296266784701563122422793621285081534865709752204207574278710231500018",
		"242.6864872978287095851298203556544621261977668634771588188722161922172718727964552121357232957194092",
		"245.8280996982398071075520620851352194267027876499343957611836025506864338750334096847784335144958054",
		"248.9697116003099371623525166327394102370078134171810019971643411931745579274551889472155318037732121",
		"252.1113230226685940010349453789353464416860540715254118029907399459605461078767505180249309866760859",
		"255.2529339830281320490805994144356008386743850488582178341457380980321641722302653315481522355301037",
		"258.3945444982395187642311505876747171215324988642865805381978593353297518140933919298534797801834920",
		"261.5361545843440692973802526071086346937871713730199924538314605607325280697133706342870162203376409",
		"264.6777642566214968097543525790789383572951299217364530760524815076745889505371920907871334845653231",
		"267.8193735296345809709430597123344754439902303862306830361991519555255862771964552049295437581134586",
		"270.9609824172707291023153605281175594834975480944740366421140058786087860466088992856241001433081492",
		"274.1025909327806792647523533889262934563161450844446786616559646749515885372734228476027571518076510",
		"277.2441990888145719904636400218613039342006275501637499909480296443081242040382185204041150697874797",
		"280.3858068974555970383791802819262494134271357718078057712268994287169495986725526885715675978294858",
		"283.5274143702514032587151817118603398632104226631419944974490323571727795725035487516958664929718895",
		"286.6690215182434431627341479961726127789883817383946477134027249375340726098980767623191531400498132",
		"289.8106283519944089128675606062573148003855915913275338996461034958131915923342775072057639746841452",
		"292.9522348816139030037286541183898679860357475415608291295070353038693289546266237925044244830760758",
		"296.0938411167824747437305359223346272879946671498746210155027390522424445245003443920186187511291827",
		"299.2354470667741426352574091585558040026118013067531043337574401524685212378464006679826400966038039",
		"302.3770527404775127692335594013707809100578831419259083403583911029237174188103967939153484710183268",
		"305.5186581464155942916187946939532935146396365807115756988835093331947269516548569472211532280304987",
		"308.6602632927644047708067635287482915627607586160043851510832939613429005356256191645569858775168825",
		"311.8018681873704508125112416665551540310176002572047531633813161412176840920137410031972570629115877",
		"314.9434728377671624580656002456121124436876018887978599031091248782240473746451805583714881611297877"
	};

	// Test with 16, 32, 64 digits precision
	int prec = 16;
	for (int i = 0; i < 3; ++i)
	{
		mpfr::mpreal::set_default_prec(mpfr::digits2bits(prec));

		std::vector<mpfr::mpreal> zeros(101);
		casema::besselZerosJ1(zeros.size(), zeros.data());

		for (int j = 0; j < ref.size(); ++j)
		{
			CAPTURE(prec);
			CAPTURE(j);
			CHECK(equal(ref[j], zeros[j]));
		}

		prec *= 2;
	}
}
