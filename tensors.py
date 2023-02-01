import numpy as np
import pandas as pd
import os
import argparse
import logging
import cooler
from collections import namedtuple
import torch
from expecto import Beluga
import pyfasta
import time

from configs import load_genes
from args import add_args, setup
from filenames import get_pred_cache_filename, get_tensor_filename

Region = namedtuple("Region", ["chrom", "start", "end"])


def encode_seqs(seqs, inputsize=2000):
    enc_seqs = np.zeros((len(seqs), 4, inputsize), np.bool_)

    enc_dict = {'A': np.asarray([1, 0, 0, 0]), 'G': np.asarray([0, 1, 0, 0]),
                'C': np.asarray([0, 0, 1, 0]), 'T': np.asarray([0, 0, 0, 1]),
                'N': np.asarray([0, 0, 0, 0]), 'H': np.asarray([0, 0, 0, 0]),
                'a': np.asarray([1, 0, 0, 0]), 'g': np.asarray([0, 1, 0, 0]),
                'c': np.asarray([0, 0, 1, 0]), 't': np.asarray([0, 0, 0, 1]),
                'n': np.asarray([0, 0, 0, 0]), '-': np.asarray([0, 0, 0, 0])}

    for n, seq in enumerate(seqs):
        for i, c in enumerate(seq):
            enc_seqs[n, :, i] = enc_dict[c]
    return enc_seqs


def predict(seqs_to_predict, window_size, batch_size, cuda, model):
    predictions = list()
    enc_seqs = encode_seqs(seqs_to_predict, inputsize=window_size).astype(np.float32)
    for e in range(0, enc_seqs.shape[0], batch_size):
        model_input = torch.from_numpy(np.array(enc_seqs[e:e + batch_size])).unsqueeze(2)
        if cuda:
            model_input = model_input.cuda()
        prediction = model.forward(model_input).cpu().detach().numpy().copy()
        predictions.append(prediction)
    return np.vstack(predictions)


def cached_region(pred_cache_folder, cached_predictions, reg):
    if reg in cached_predictions:
        return cached_predictions[reg]
    else:
        cached_pred_filename = get_pred_cache_filename(pred_cache_folder, reg.chrom, reg.start, reg.end)
        if os.path.isfile(cached_pred_filename):
            try:
                tensor = np.load(cached_pred_filename)["pred"]
                cached_predictions[reg] = tensor
                return tensor
            except ValueError:
                pass
    return None


def predict_regions_with_cache(pred_cache_folder, genome, model, cached_predictions, regs_to_predict, window_size,
                               features_count, batch_size, cuda):
    start_time = time.time()
    tensor = np.ndarray(shape=(len(regs_to_predict), features_count), dtype=np.float32)
    cache_used = 0
    seqs_to_predict = dict()
    for i, reg in enumerate(regs_to_predict):
        cached_tensor = cached_region(pred_cache_folder, cached_predictions, reg)
        if cached_tensor is not None:
            tensor[i] = cached_tensor
            cache_used += 1
        else:
            seq = genome.sequence({'chr': reg.chrom, 'start': reg.start, 'stop':  reg.end})
            seqs_to_predict[i] = seq
    if len(seqs_to_predict):
        predictions = predict(seqs_to_predict.values(), window_size, batch_size, cuda, model)
        for p, i in enumerate(seqs_to_predict.keys()):
            reg = regs_to_predict[i]
            cached_predictions[reg] = predictions[p]
            tensor[i] = predictions[p]
            cached_pred_filename = get_pred_cache_filename(pred_cache_folder, reg.chrom, reg.start, reg.end)
            np.savez_compressed(cached_pred_filename, pred=predictions[p])
    if len(regs_to_predict) > 0:
        logging.info(f"{len(regs_to_predict)/(time.time() - start_time):.2f} bins/s. Cache used: {cache_used}"
                     f" ({100*cache_used/len(regs_to_predict):.2f}%)")
    return tensor

def naive_regions_and_pos_weights(gene, window_size, bin_size, span=20000): 
    regs_to_predict, pos_weights = expecto_regions_and_pos_weights(gene, window_size, bin_size, span)
    pos_weights[True] = 1
    pos_weights = np.vstack([pos_weights[0], pos_weights[0]]) # we take only 2 positions to ensure same size of tensor as in case Hi-C was present; weights are 1 anyway so it doesn't matter what we take
    return regs_to_predict, pos_weights


def expecto_regions_and_pos_weights(gene, window_size, bin_size, span=20000):
    regs_to_predict = list()
    pos_weight_shifts = np.array(list(range(-span, span, bin_size)))
    for shift in pos_weight_shifts:
        chrom, tss, strand = gene.seqnames, gene.CAGE_representative_TSS, 1 if gene.strand == "+" else -1
        reg = Region(chrom, tss + shift * strand - int(0.5*window_size - 1), tss + shift * strand + int(0.5*window_size))
        regs_to_predict.append(reg)

    pos_weights = np.vstack([
        np.exp(-0.01*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.02*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.05*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.1*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.2*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts <= 0),
        np.exp(-0.01*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0),
        np.exp(-0.02*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0),
        np.exp(-0.05*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0),
        np.exp(-0.1*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0),
        np.exp(-0.2*np.abs(pos_weight_shifts)/200)*(pos_weight_shifts >= 0)])

    return regs_to_predict, pos_weights


def hic_regions_and_pos_weights(gene, mcools, hic_resolution, max_dist, bin_size, window_size, tensor_width, progress):
    chrom = gene.seqnames
    strand = gene.strand
    tss = hic_resolution * round(gene.CAGE_representative_TSS / hic_resolution)
    assert(all(mcools[0].offset(f"{chrom}:{tss}-{tss}") == mcool.offset(f"{chrom}:{tss}-{tss}") for mcool in mcools))

    tss_bin = mcools[0].offset(f"{chrom}:{tss}-{tss}")
    max_span = 10**10
    start, end = max(tss - max_span, 0), min(tss + max_span, mcools[0].chromsizes[chrom])
    counts = list()

    for mcool in mcools:
        counts.append(mcool.matrix(field="count", balance=None, as_pixels=True, ignore_index=True, join=False)
                      .fetch(f"{chrom}:{start}-{tss}", f"{chrom}:{tss}-{tss+hic_resolution}"))
        counts.append(mcool.matrix(field="count", balance=None, as_pixels=True, ignore_index=True, join=False)
                      .fetch(f"{chrom}:{tss}-{tss+hic_resolution}", f"{chrom}:{tss}-{end}"))
    counts = pd.concat(counts)
    counts = counts.groupby(["bin1_id", "bin2_id"])['count'].sum().reset_index()

    logging.info(f"{gene.id} {gene.symbol} {gene.seqnames}:{gene.CAGE_representative_TSS} cnts: {len(counts)} {progress:.2f}%")

    if len(counts) == 0:
        logging.info(f"{gene.id} Not enough 3D data: {len(counts)} {progress:.2f}%")
        return None, None

    # determine seq positions and transformation
    tensor_width2 = tensor_width // 2

    pos_dists = list()
    for count in counts.itertuples():
        dist = 1 / (count.count**0.5)
        if dist >= max_dist:
            continue
        x_bin = count.bin2_id if count.bin1_id == tss_bin else count.bin1_id
        pos = tss + hic_resolution * (x_bin - tss_bin)
        if abs(pos - tss) < 20000:
            continue
        pos_dists.append((pos, dist))

    regs_to_predict, pos_weights = list(), np.zeros(shape=(tensor_width, (hic_resolution//bin_size) * len(pos_dists)))

    i = 0
    for pos, dist in pos_dists:
        for j in range(hic_resolution//bin_size):
            pos_bin = pos + j * bin_size
            reg = Region(chrom, pos - int(0.5*window_size - 1), pos + int(0.5*window_size))
            regs_to_predict.append(reg)
            second_half = (strand == "+" and pos_bin < tss) or (strand == "-" and pos_bin >= tss)
            pos_weights[tensor_width2 * second_half, i] = 1
            i += 1

    return regs_to_predict, pos_weights


def generate_tensors(args):
    # load protein coding genes and Chromatin Contact Domains (CCDs), find CCDs containing genes
    genes = load_genes(args.genes_filename, args.gene_type, args.data_folder)

    mcools = None
    if args.mcool_filename:
        mcool_filename = os.path.expanduser(args.mcool_filename)
        mcool_hic = cooler.Cooler(f"{mcool_filename}::/resolutions/{args.hic_resolution}")
        mcools = [mcool_hic]

    logging.info(f"Loading genome from {args.genome_filename}")
    genome = pyfasta.Fasta(args.genome_filename)

    logging.info(f"Done. Loading Beluga model from {args.beluga_model_filename}")
    model = Beluga()
    model.load_state_dict(torch.load(args.beluga_model_filename))
    model.eval()
    if args.compute_platform == "CUDA":
        model.cuda()
    logging.info("Done.")

    if args.divisor is not None and args.reminder is not None:
        genes = genes[(genes.index % args.divisor) == args.reminder]
    if args.test_only:
        genes = genes[genes.seqname == "chr8"]
    cached_preds = dict()

    # for each gene, calculate distances from TSSs and build tensor
    for i, gene in enumerate(genes.itertuples()):
        logging.info(f"{gene.id} {gene.symbol} {gene.seqnames}:{gene.CAGE_representative_TSS} {i+1}/{len(genes)}")
        tensor_filename = get_tensor_filename(gene, args.data_folder)
        if not os.path.isfile(tensor_filename):
            expecto_regs, expecto_pos_weights = expecto_regions_and_pos_weights(gene, args.window_size, args.bin_size)
            expecto_tensor = predict_regions_with_cache(args.pred_cache_folder, genome, model, cached_preds, expecto_regs,
                                                        args.window_size, args.features_count, args.batch_size,
                                                        args.compute_platform == "CUDA")
            expecto_tensor = np.sum(expecto_pos_weights[:, :, None]*expecto_tensor[None, :, :], axis=1)
            expecto_tensor = expecto_tensor.flatten()

            if mcools:
                hic_regions, hic_pos_weights = hic_regions_and_pos_weights(gene, mcools, args.hic_resolution, args.max_dist,
                                                                           args.bin_size, args.window_size, 2,
                                                                           100*i/len(genes))
                if hic_regions is None:
                    logging.info("No hic_regions")
                    hic_regions, hic_pos_weights = naive_regions_and_pos_weights(gene, args.window_size, args.bin_size)
                else:
                    logging.info("Hic_regions present")
                hic_tensor = predict_regions_with_cache(args.pred_cache_folder, genome, model, cached_preds, hic_regions,
                                                        args.window_size, args.features_count, args.batch_size, 
                                                        args.compute_platform == "CUDA")
                hic_tensor = np.sum(hic_pos_weights[:, :, None]*hic_tensor[None, :, :], axis=1)
                hic_tensor = hic_tensor.flatten()
                tensor = np.concatenate([expecto_tensor, hic_tensor])
            else:
                tensor = expecto_tensor

            np.save(tensor_filename, tensor)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_args(parser)
    args = parser.parse_args()
    setup(args)
    generate_tensors(args)
