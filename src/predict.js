'use strict';

function testOnData(sequence, threshold) {

    /*run on data*/
    var result = [];
    sequence.forEach(function(arrayItem) {
        var output = net.run(arrayItem);
        // console.log(output);
        result.push(output);
    });

    /*Calc % of true positives using a threshold
    threshold must be between 0 and 1*/
    const filterAlpha = result.filter(property => property.alpha >= threshold);
    const filterBeta = result.filter(property => property.beta >= threshold);

    var trueAlpha = filterAlpha.length / result.length;
    var roundedTrueAlpha = Math.round(trueAlpha * 100) / 100;
    var trueBeta = filterBeta.length / result.length;
    var roundedTrueBeta = Math.round(trueBeta * 100) / 100;

    console.log(`Test pour alpha`);
    console.log(`True positives : ${roundedTrueAlpha * 100}%`);
    console.log(`Test pour beta`);
    console.log(`True positives : ${roundedTrueBeta * 100}%`);

}

/*Find features and store them*/

function process(data) {
    let result = [];

    /* Read data from fasta format*/
    let tab = data.split("\n");

    for (let i = 0; i < tab.length; i++) {
        let sequence = "";
        if (tab[i] != "" && tab[i][0] != ">") { //remove "\n" and concatenate untill reaches a new ">" 
            while (tab[i] != "" && tab[i][0] != ">") {
                sequence += tab[i];
                i++;
            }
            result.push(sequence);
        }
    }

    /*Hydrophobicity Kyte-Doolittle scale
    Isoleucine 	I 	4.5
    Valine 	V 	4.2
    Leucine 	L 	3.8
    Phenylalanine 	F 	2.8
    Cysteine 	C 	2.5
    Methionine 	M 	1.9
    Alanine 	A 	1.8
    Glycine 	G 	-0.4
    Threonine 	T 	-0.7
    Serine 	S 	-0.8
    Tryptophan 	W 	-0.9
    Tyrosine 	Y 	-1.3
    Proline 	P 	-1.6
    Histidine 	H 	-3.2
    Glutamic acid 	E 	-3.5
    Glutamine 	Q 	-3.5
    Aspartic acid 	D 	-3.5
    Asparagine 	N 	-3.5
    Lysine 	K 	-3.9
    Arginine 	R 	-4.5
    */

    let KDHydrophobicity = { "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5, "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5 };
    let H = [];
    let averageH = 0;

    for (let i = 0; i < result.length; i++) {
        let hydropathie = [];
        for (let j = 0; j < result[i].length; j++) {
            if (j > 2 && j < result[i].length - 3) {
                averageH = (KDHydrophobicity[result[i][j - 3]] + KDHydrophobicity[result[i][j - 2]] + KDHydrophobicity[result[i][j - 1]] +
                    KDHydrophobicity[result[i][j]] + KDHydrophobicity[result[i][j + 1]] + KDHydrophobicity[result[i][j + 2]] + KDHydrophobicity[result[i][j + 3]]) / 7;
                hydropathie.push(averageH);
            } else {
                averageH = KDHydrophobicity[result[i][j]];
                hydropathie.push(averageH);
            }
        }
        let averageHydropathie = 0;
        for (let k = 0; k < hydropathie.length; k++) {
            averageHydropathie += hydropathie[k];
        }
        averageHydropathie /= hydropathie.length;
        H.push(averageHydropathie)

    }

    /*Calc features for the sequences in the fasta file
    Infos : AA impliqués dans la formation des hélices alpha : MALEK
            AA impliqués dans la formation des feuillets beta: FYWTVIL et P aux extrémités.
    */
    let features = [];
    for (let i = 0; i < result.length; i++) {
        let matrix = { "A": 0, "D": 0, "R": 0, "N": 0, "C": 0, "E": 0, "Q": 0, "G": 0, "H": 0, "I": 0, "L": 0, "K": 0, "M": 0, "F": 0, "P": 0, "S": 0, "T": 0, "W": 0, "Y": 0, "V": 0 };
        for (let j = 0; j < result[i].length; j++) {
            matrix[result[i][j]] += 1;
        }
        let len = result[i].length;

        let aromaticity = (matrix["H"] + matrix["F"] + matrix["W"] + matrix["Y"]) / len;
        matrix["aromaticity"] = aromaticity;

        for (let key in matrix) {
            matrix[key] = matrix[key] / len;
        }
        matrix["hydrophibicity"] = H[i];


        features[i] = matrix;
    }
    return features;

}


/*MAIN*/
var betaFeatures = process(dataset);
var alphaTest = process(alpha);
var alphaHTest = process(alphaH);

var test2 = process(Test2);
var test3 = process(Test3);
var test4 = process(Test4);
var test5 = process(Test5);


// /*Training and inferences*/
// /* Features format : [{input:{features} , output:{beta:1,alpha:0}}]*/
var net = new brain.NeuralNetwork();
var config = [];

/*Train classifier on alphaG & beta*/
test2.forEach(function(arrayItem) {
    config.push({ input: arrayItem, output: { beta: 1, alpha: 0 } });
});

test3.forEach(function(arrayItem) {
    config.push({ input: arrayItem, output: { beta: 1, alpha: 0 } });
});


alphaTest.forEach(function(arrayItem) {
    config.push({ input: arrayItem, output: { beta: 0, alpha: 1 } });
});

alphaHTest.forEach(function(arrayItem) {
    config.push({ input: arrayItem, output: { beta: 0, alpha: 1 } });
});

console.log("Training has begun ! Wait...");
net.train(config, {
    // Defaults values --> expected validation
    iterations: 20000, // the maximum times to iterate the training data --> number greater than 0
    errorThresh: 0.005, // the acceptable error percentage from training data --> number between 0 and 1
    log: true, // true to use console.log, when a function is supplied it is used --> Either true or a function
    logPeriod: 10, // iterations between logging out --> number greater than 0
    learningRate: 0.3, // scales with delta to effect training rate --> number between 0 and 1
    momentum: 0.1, // scales with next layer's change value --> number between 0 and 1
    callback: null, // a periodic call back that can be triggered while training --> null or function
    callbackPeriod: 10, // the number of iterations through the training data between callback calls --> number greater than 0
    timeout: Infinity // the max number of milliseconds to train for --> number greater than 0
});
// net.train(config);
// console.log(net.train(config));
console.log("Training, done. Yay!");

// /*Test classifier*/
let threshold = 0.9;

testOnData(betaFeatures, threshold);
console.log("-------------------------------------------------------------------");
testOnData(test4, threshold);
console.log("-------------------------------------------------------------------");
testOnData(test5, threshold);