"""Exactly 30 seconds of original 140-BPM score and concise Chinese narration."""
from pathlib import Path
import argparse
import asyncio
import hashlib
import json
import os
import shutil
import subprocess
import wave
import edge_tts
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RATE = 48000


def tool(name):
    found = shutil.which(name)
    if not found:
        raise FileNotFoundError(name)
    return str(Path(found).resolve())


def run(command):
    result = subprocess.run(command, capture_output=True)
    if result.returncode:
        raise RuntimeError(result.stderr.decode('utf-8', errors='replace')[-1500:])
    return result.stdout


def decode(path, speed=1):
    data = run([tool('ffmpeg'), '-v', 'error', '-i', str(path), '-af', f'atempo={speed:.8f}',
                '-f', 'f32le', '-ac', '1', '-ar', str(RATE), 'pipe:1'])
    return np.frombuffer(data, dtype='<f4').copy()


def cache_path(scene, story, work):
    identity = '\n'.join((story['voice'], story['voice_rate'], scene['narration']))
    key = hashlib.sha256(identity.encode('utf-8')).hexdigest()[:12]
    return work / f"voice-{scene['id']}-{key}.mp3"


async def speech(scene, story, work, semaphore):
    if not scene['narration']:
        return
    path = cache_path(scene, story, work)
    if path.exists() and path.stat().st_size > 1000:
        return
    async with semaphore:
        await edge_tts.Communicate(scene['narration'], story['voice'], rate=story['voice_rate']).save(str(path))
        print('Voice ready: ' + scene['id'], flush=True)


async def voices(story, work):
    semaphore = asyncio.Semaphore(2)
    await asyncio.gather(*(speech(scene, story, work, semaphore) for scene in story['scenes']))


def score(duration, bpm, cuts):
    """Procedural synth: original notes, percussion and sound design, no samples."""
    rng = np.random.default_rng(140307)
    count = round(duration * RATE)
    output = np.zeros((count, 2), dtype=np.float32)
    beat = 60 / bpm

    def add(signal, at, gain, pan=0):
        start = max(0, round(at * RATE))
        length = min(len(signal), count - start)
        if length <= 0:
            return
        output[start:start + length, 0] += signal[:length] * gain * (1 - pan * .35)
        output[start:start + length, 1] += signal[:length] * gain * (1 + pan * .35)

    def time(length):
        return np.arange(round(length * RATE), dtype=np.float32) / RATE

    chords = [(38, 50, 53, 57), (34, 46, 50, 53), (41, 53, 57, 60), (36, 48, 52, 55)]
    for index in range(70):
        at = index * beat
        chord = chords[(index // 8) % len(chords)]
        t = time(.24)
        phase = 2 * np.pi * (52 * t + 3.2 * (1 - np.exp(-38 * t)))
        add(np.sin(phase) * np.exp(-21 * t), at, .15)
        if index % 2:
            t = time(.12)
            noise = rng.standard_normal(len(t)).astype(np.float32)
            add((noise - np.convolve(noise, np.ones(7) / 7, mode='same')) * np.exp(-35 * t), at, .035)
        t = time(beat * .82)
        hz = 440 * 2 ** ((chord[0] - 69) / 12)
        bass = (np.sin(2 * np.pi * hz * t) + .24 * np.sin(4 * np.pi * hz * t)) * np.sin(np.pi * t / (beat * .82))
        add(bass, at + .018, .075)
        for sub in range(4):
            t = time(.055)
            noise = rng.standard_normal(len(t)).astype(np.float32)
            high = noise - np.convolve(noise, np.ones(11) / 11, mode='same')
            add(high * np.exp(-95 * t), at + sub * beat / 4, .013 if sub % 2 else .02, (-1) ** sub * .7)
        for sub in range(2):
            t = time(.19)
            midi = chord[1 + (index + sub) % 3] + (12 if index % 4 == 3 else 0)
            hz = 440 * 2 ** ((midi - 69) / 12)
            synth = np.sin(2 * np.pi * hz * t) + .2 * np.sin(2 * np.pi * hz * 2.003 * t)
            add(synth * np.exp(-22 * t), at + sub * beat / 2, .04, (-1) ** index * .55)
    for at in cuts:
        t = time(.18)
        noise = rng.standard_normal(len(t)).astype(np.float32)
        add(noise * np.sin(np.pi * t / .18) ** 2, max(0, at - .13), .032, -.3)
        t = time(.2)
        add(np.sin(2 * np.pi * (87 * t + 1.5 * (1 - np.exp(-30 * t)))) * np.exp(-23 * t), at, .08)
    t = np.arange(count, dtype=np.float32) / RATE
    fade = np.minimum(t / .08, 1) * np.minimum(np.maximum(duration - t, 0) / .5, 1)
    return output * fade[:, None]


def stamp(value):
    ms = round(value * 1000)
    return f'{ms // 3600000:02}:{ms // 60000 % 60:02}:{ms // 1000 % 60:02}.{ms % 1000:03}'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--direct', action='store_true')
    parser.add_argument('--work-dir', type=Path, default=ROOT / '.local/video-production/v2/audio')
    args = parser.parse_args()
    work = args.work_dir.resolve()
    if not work.is_relative_to(ROOT):
        raise ValueError('Keep scratch inside the repository')
    work.mkdir(parents=True, exist_ok=True)
    for key in ('TEMP', 'TMP', 'TMPDIR'):
        os.environ[key] = str(work)
    if args.direct:
        for key in list(os.environ):
            if key.lower() in ('http_proxy', 'https_proxy', 'all_proxy'):
                os.environ.pop(key)
        os.environ['NO_PROXY'] = os.environ['no_proxy'] = '*'
    story = json.loads((HERE / 'storyboard.json').read_text(encoding='utf-8'))
    asyncio.run(voices(story, work))
    narration = np.zeros(round(story['duration'] * RATE), dtype=np.float32)
    timing, cues = [], []
    for scene in story['scenes']:
        start = scene['start_frame'] / story['fps']
        end = scene['end_frame'] / story['fps']
        row = {**scene, 'start': round(start, 6), 'duration': round(end - start, 6)}
        if scene['narration']:
            path = cache_path(scene, story, work)
            voice = decode(path)
            available = scene.get('voice_end_frame', scene['end_frame']) / story['fps'] - start - .2
            factor = max(1, (len(voice) / RATE) / available)
            if factor > 1.4:
                raise ValueError(f"Shorten {scene['id']} narration; compression would be {factor:.2f}x")
            if factor > 1:
                voice = decode(path, factor)
            offset = round((start + .12) * RATE)
            narration[offset:offset + len(voice)] += voice
            row.update(voice_start=round(offset / RATE, 6), voice_duration=round(len(voice) / RATE, 6), tempo_factor=round(factor, 6))
            cues.append((offset / RATE, (offset + len(voice)) / RATE, scene['narration'].replace('M六', 'M6')))
        timing.append(row)
    soundtrack = score(story['duration'], story['bpm'], [s['start'] for s in timing[1:]])
    absolute = np.abs(narration).astype(np.float64)
    integral = np.concatenate(([0.], np.cumsum(absolute)))
    indices = np.arange(len(narration))
    beginnings = np.maximum(indices - round(.1 * RATE), 0)
    mean = (integral[indices + 1] - integral[beginnings]) / (indices + 1 - beginnings)
    duck = np.where(mean > .01, .36, .85).astype(np.float32)
    mix = narration[:, None] * .95 + soundtrack * duck[:, None]
    peak = float(np.abs(mix).max())
    if peak > .95:
        mix *= .95 / peak
    wav = work / 'master-mix.wav'
    with wave.open(str(wav), 'wb') as output:
        output.setnchannels(2)
        output.setsampwidth(2)
        output.setframerate(RATE)
        output.writeframes((mix * 32767).astype('<i2').tobytes())
    master = HERE / 'assets/master.m4a'
    run([tool('ffmpeg'), '-v', 'error', '-y', '-i', str(wav), '-af', 'loudnorm=I=-14:TP=-1.5:LRA=7',
         '-ar', str(RATE), '-c:a', 'aac', '-b:a', '160k', '-t', '30', '-movflags', '+faststart', str(master)])
    data = {key: story[key] for key in ('duration', 'fps', 'bpm', 'voice', 'voice_rate')}
    data.update(schema_version=2, scenes=timing)
    (HERE / 'timing.json').write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    caption = ROOT / 'assets/video/engineering-research-prototypes-intro.zh.vtt'
    caption.write_text('WEBVTT\n\n' + '\n\n'.join(f'{i + 1}\n{stamp(a)} --> {stamp(b)}\n{text}' for i, (a, b, text) in enumerate(cues)) + '\n', encoding='utf-8', newline='\n')
    print(json.dumps({'seconds': story['duration'], 'cues': len(cues), 'scenes': len(timing), 'voice_fits': True}), flush=True)


if __name__ == '__main__':
    main()
