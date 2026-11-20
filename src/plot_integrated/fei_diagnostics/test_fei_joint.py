"""Check boundaries and heterogeneous corrections in the joint diagnostics."""
import numpy as np
import pandas as pd
import pytest

from analyze_fei_joint import application_table, joint_table


def test_probability_boundary_and_missing_truth_remain_distinct():
    events = pd.DataFrame([
        dict(category=3, sig_mode=0, sample='DATA', weight=1., tag_probability=.01, tag_is_signal=np.nan),
        dict(category=3, sig_mode=0, sample='SIGNAL', weight=2., tag_probability=.01, tag_is_signal=1.),
        dict(category=3, sig_mode=0, sample='SIGNAL', weight=3., tag_probability=.02, tag_is_signal=np.nan),
        dict(category=3, sig_mode=0, sample='SIGNAL', weight=4., tag_probability=np.nan, tag_is_signal=0.),
    ])
    table = joint_table(events).set_index('probability_region')
    assert table.loc['le_0.01', 'data'] == 1
    assert table.loc['le_0.01', 'mc'] == 2
    assert table.loc['gt_0.01', 'data'] == 0
    assert table.loc['gt_0.01', 'mc_unknown'] == 3
    assert table.loc['gt_0.01', 'ratio_error_high'] > 0
    assert table.loc['unknown', 'mc_incorrect'] == 4
    assert table.mc_sumw2.sum() == 29


def test_correct_only_application_keeps_official_weights_on_unknown():
    # Unequal weights and official factors expose an accidental global scaling
    # or accidental treatment of unknown truth as correctly reconstructed.
    events = pd.DataFrame([
        dict(category=3, sample='DATA', weight=1., tag_is_signal=np.nan, official_fei=1.),
        dict(category=3, sample='SIGNAL', weight=2., tag_is_signal=1., official_fei=.8),
        dict(category=3, sample='SIGNAL', weight=3., tag_is_signal=0., official_fei=.5),
        dict(category=3, sample='SIGNAL', weight=5., tag_is_signal=np.nan, official_fei=.2),
    ])
    result = application_table(events).iloc[0]
    assert result.inclusive_ratio == pytest.approx(.1)
    assert result.mc_after_ratio_on_all == pytest.approx(1.)
    assert result.official_mc_noncorrect_including_unknown == pytest.approx(2.5)
    assert result.mc_after_ratio_on_correct_official_on_rest == pytest.approx(2.7)
