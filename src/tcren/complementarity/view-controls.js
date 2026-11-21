// Shared controls for static molecular and table views.
const nav = document.getElementById('view-controls');
const theme = document.createElement('select');
for (const [value, label] of [['system','System'], ['light','Light'], ['dark','Dark']]) {
    theme.add(new Option(label, value));
}
const label = document.createElement('label');
label.append('Theme ', theme);
nav.append(label);
const preference = matchMedia('(prefers-color-scheme: dark)');
theme.value = localStorage.getItem('tcren-dashboard-theme') || 'system';
function applyTheme() {
    document.documentElement.dataset.theme = theme.value === 'system'
        ? (preference.matches ? 'dark' : 'light') : theme.value;
}
applyTheme();
preference.addEventListener('change', applyTheme);
theme.addEventListener('change', () => {
    localStorage.setItem('tcren-dashboard-theme', theme.value);
    applyTheme();
});

const response = await fetch('navigation.json').catch(() => null);
if (response?.ok) {
    const {links = [], corpora = []} = await response.json();
    if (corpora.length) {
        const corpus = document.createElement('select');
        for (const item of corpora) {
            const option = new Option(item.label, item.href, item.current, item.current);
            corpus.add(option);
        }
        corpus.addEventListener('change', () => location.assign(corpus.value));
        const chooser = document.createElement('label');
        chooser.append('Corpus ', corpus);
        nav.prepend(chooser);
    }
    for (const item of links) {
        const link = document.createElement('a');
        link.textContent = item.label;
        if (item.current) link.setAttribute('aria-current', 'page');
        const sync = () => {
            const target = new URL(item.href, location.href);
            if (item.model) {
                const current = new URL(location.href);
                for (const key of ['group', 'id']) {
                    if (current.searchParams.has(key)) target.searchParams.set(key, current.searchParams.get(key));
                }
            }
            link.href = target.href;
        };
        sync();
        document.addEventListener('tcren:model', sync);
        nav.append(link);
    }
}
